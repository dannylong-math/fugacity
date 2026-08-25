import { access, readdir, readFile } from "node:fs/promises";
import { dirname, extname, relative, resolve, sep } from "node:path";
import { fileURLToPath } from "node:url";

const docsDirectory = dirname(fileURLToPath(import.meta.url));
const outputDirectory = resolve(docsDirectory, "dist");
const xmlDirectory = resolve(docsDirectory, "../build/doxygen/xml");

async function filesBelow(directory) {
  const files = [];
  for (const entry of await readdir(directory, { withFileTypes: true })) {
    const path = resolve(directory, entry.name);
    if (entry.isDirectory()) files.push(...await filesBelow(path));
    else if (entry.isFile()) files.push(path);
  }
  return files;
}

function insideOutput(path) {
  return path === outputDirectory || path.startsWith(`${outputDirectory}${sep}`);
}

function idsIn(html) {
  return new Set([...html.matchAll(/\bid=["']([^"']+)["']/g)].map((match) => match[1]));
}

const allFiles = await filesBelow(outputDirectory);
const artifacts = new Set(allFiles);
const htmlFiles = allFiles.filter((path) => extname(path) === ".html");
const htmlCache = new Map();
let checkedLinks = 0;
const failures = [];

for (const source of htmlFiles) {
  const html = await readFile(source, "utf8");
  htmlCache.set(source, html);
  for (const match of html.matchAll(/\b(?:href|src)=["']([^"']+)["']/gi)) {
    const rawHref = match[1];
    if (/^(?:[a-z][a-z0-9+.-]*:|\/\/)/i.test(rawHref)) continue;
    const [pathAndQuery, fragment = ""] = rawHref.split("#", 2);
    const pathPart = decodeURIComponent(pathAndQuery.split("?", 1)[0]);
    let target;
    if (pathPart === "") target = source;
    else if (pathPart.startsWith("/")) target = resolve(outputDirectory, `.${pathPart}`);
    else target = resolve(dirname(source), pathPart);
    if (pathPart.endsWith("/")) target = resolve(target, "index.html");
    checkedLinks += 1;
    if (!insideOutput(target) || !artifacts.has(target)) {
      failures.push(`${relative(outputDirectory, source)} -> ${rawHref}`);
      continue;
    }
    if (fragment && extname(target) === ".html") {
      const targetHtml = htmlCache.get(target) ?? await readFile(target, "utf8");
      htmlCache.set(target, targetHtml);
      if (!idsIn(targetHtml).has(decodeURIComponent(fragment))) {
        failures.push(`${relative(outputDirectory, source)} -> ${rawHref} (missing fragment)`);
      }
    }
  }
}

if (failures.length > 0) {
  throw new Error(
    `${failures.length} exact static artifact links failed:\n${failures.slice(0, 40).join("\n")}`,
  );
}

for (const legacy of [
  "getting_started.html",
  "tutorial.html",
  "implementing_a_new_eos.html",
  "api/concepts.html",
]) {
  await access(resolve(outputDirectory, legacy));
}

const coreXml = await readFile(resolve(xmlDirectory, "group__core.xml"), "utf8");
const publicMembers = [...coreXml.matchAll(
  /<memberdef\b[^>]*\bkind="(function|variable)"[^>]*>[\s\S]*?<name>([^<]+)<\/name>[\s\S]*?<\/memberdef>/g,
)];
const functions = publicMembers.filter((match) => match[1] === "function");
const variables = publicMembers.filter((match) => match[1] === "variable");
if (functions.length === 0 || !variables.some((match) => match[2] === "ideal_gas_constant")) {
  throw new Error("Core API inventory lacks public free functions or ideal_gas_constant.");
}
const corePage = resolve(outputDirectory, "api/core.html");
const coreHtml = await readFile(corePage, "utf8");
const renderedMemberHeadings = [...coreHtml.matchAll(/<h3\b[^>]*\bid=["'][^"']+["']/g)];
if (renderedMemberHeadings.length !== publicMembers.length) {
  throw new Error(
    `Core API rendered ${renderedMemberHeadings.length} member sections for `
      + `${publicMembers.length} public function/variable definitions.`,
  );
}
for (const name of new Set(publicMembers.map((match) => match[2]))) {
  if (!coreHtml.includes(`>${name}</h3>`)) {
    throw new Error(`Public core member ${name} is absent from api/core.html.`);
  }
}

for (const concept of ["EquationOfState", "IdealEoS", "ResidualEoS"]) {
  const conceptPath = resolve(outputDirectory, `api/fugacity-${concept}.html`);
  const conceptHtml = await readFile(conceptPath, "utf8");
  if (!conceptHtml.includes(concept)) throw new Error(`Concept ${concept} is absent.`);
}

const apiHtml = [...htmlCache.entries()]
  .filter(([path]) => relative(outputDirectory, path).startsWith(`api${sep}`))
  .map(([, html]) => html)
  .join("\n");
if (/fugacity::detail|__enzyme|Enzyme/i.test(apiHtml)) {
  throw new Error("Internal fugacity::detail or Enzyme declarations leaked into the site.");
}
const combinedHtml = [...htmlCache.values()].join("\n");
if (/Edit this page/.test(combinedHtml)) {
  throw new Error("Generated site contains an unverified Edit-this-page link.");
}

process.stdout.write(
  `Checked ${checkedLinks} exact static artifact links across ${htmlFiles.length} HTML pages; `
    + `verified ${functions.length} core function overloads, ideal_gas_constant, `
    + "3 concepts, 4 legacy URLs, and no internal API leakage.\n",
);
