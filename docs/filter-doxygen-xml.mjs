import { readdir, readFile, writeFile } from "node:fs/promises";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const docsDirectory = dirname(fileURLToPath(import.meta.url));
const xmlDirectory = resolve(process.argv[2] ?? resolve(docsDirectory, "../build/doxygen/xml"));
const xmlFiles = (await readdir(xmlDirectory))
  .filter((name) => name.endsWith(".xml"));

const privateMemberIds = new Set();
const privateMember = /\s*<memberdef\b[^>]*\bprot="private"[^>]*>[\s\S]*?<\/memberdef>/g;
const formula = /<formula\b[^>]*>([\s\S]*?)<\/formula>/g;
const inheritance = /<(basecompoundref|derivedcompoundref)(\b[^>]*)>([\s\S]*?)<\/\1>/g;
const publicConcepts = [
  ["conceptfugacity_1_1EquationOfState", "fugacity::EquationOfState"],
  ["conceptfugacity_1_1IdealEoS", "fugacity::IdealEoS"],
  ["conceptfugacity_1_1ResidualEoS", "fugacity::ResidualEoS"],
];

function encodeFormula(entry, contents) {
  const trimmed = contents.trim();
  const display = trimmed.startsWith("\\[") && trimmed.endsWith("\\]");
  const inline = trimmed.startsWith("$") && trimmed.endsWith("$");
  if (!display && !inline) {
    throw new Error(`Unsupported Doxygen formula delimiters: ${trimmed}`);
  }

  const tex = display
    ? trimmed.slice(2, -2).trim()
    : trimmed.slice(1, -1).trim();
  const payload = Buffer.from(tex, "utf8").toString("base64url");
  return `FUGACITYKATEX:${display ? "DISPLAY" : "INLINE"}:${payload}:FUGACITYEND`;
}

function simplifyInheritance(entry, element, attributes, contents) {
  const typeName = contents.split("&lt;", 1)[0];
  return `<${element}${attributes}>${typeName}</${element}>`;
}

for (const name of xmlFiles) {
  const path = resolve(xmlDirectory, name);
  const xml = await readFile(path, "utf8");
  let filtered = xml
    .replace(privateMember, (definition) => {
      const id = definition.match(/\bid="([^"]+)"/)?.[1];
      if (id) privateMemberIds.add(id);
      return "";
    })
    .replace(formula, encodeFormula)
    .replace(inheritance, simplifyInheritance);
  if (name === "group__core.xml") {
    const missingConcepts = publicConcepts.filter(([refid]) =>
      !filtered.includes(`innerclass refid="${refid}"`));
    if (missingConcepts.length > 0) {
      const references = missingConcepts.map(([refid, name]) =>
        `    <innerclass refid="${refid}" prot="public">${name}</innerclass>`,
      ).join("\n");
      filtered = filtered.replace("  </compounddef>", `${references}\n  </compounddef>`);
    }
  }
  await writeFile(path, filtered);
}

const indexPath = resolve(xmlDirectory, "index.xml");
let index = await readFile(indexPath, "utf8");
index = index.replace(/\s*<member\b[^>]*\brefid="([^"]+)"[^>]*>.*?<\/member>/g,
  (entry, refid) => privateMemberIds.has(refid) ? "" : entry);
// Moxygen assigns each member reference to its last index.xml owner. Doxygen
// lists grouped public members under both the namespace and their group, so the
// namespace duplicate otherwise steals ownership and leaves the group empty.
let duplicateGroupMembers = 0;
index = index.replace(
  /(<compound\b[^>]*\bkind="namespace"[^>]*>[\s\S]*?<\/compound>)/g,
  (namespace) => namespace.replace(
    /\s*<member\b[^>]*\brefid="group__[^"]+"[^>]*>.*?<\/member>/g,
    () => {
      duplicateGroupMembers += 1;
      return "";
    },
  ),
);
await writeFile(indexPath, index);

const definedIds = new Set();
for (const name of xmlFiles) {
  const xml = await readFile(resolve(xmlDirectory, name), "utf8");
  for (const match of xml.matchAll(/\bid="([^"]+)"/g)) {
    definedIds.add(match[1]);
  }
}

let suppressedMemberReferences = 0;
let orphanReferences = 0;
const reference = /<ref\b([^>]*\brefid="([^"]+)"[^>]*)>([\s\S]*?)<\/ref>/g;
for (const name of xmlFiles) {
  const path = resolve(xmlDirectory, name);
  const xml = await readFile(path, "utf8");
  const filtered = xml.replace(reference, (entry, attributes, refid, content) => {
    if (attributes.includes('kindref="member"')) {
      suppressedMemberReferences += 1;
      return content;
    }
    if (definedIds.has(refid)) return entry;
    orphanReferences += 1;
    return content;
  });
  await writeFile(path, filtered);
}

process.stdout.write(
  `Filtered ${privateMemberIds.size} private implementation members and `
    + `${suppressedMemberReferences} ambiguous member references from Doxygen XML; `
    + `assigned ${duplicateGroupMembers} duplicate namespace references to their `
    + `public groups and attached ${publicConcepts.length} public concepts `
    + `(${orphanReferences} orphan references).\n`,
);

if (orphanReferences > 0) {
  throw new Error(
    `Doxygen XML contains ${orphanReferences} orphan references; refusing to build.`,
  );
}
