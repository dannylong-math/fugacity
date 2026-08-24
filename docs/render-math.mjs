import { cp, mkdir, readdir, readFile, writeFile } from "node:fs/promises";
import { dirname, relative, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import katex from "katex";

const docsDirectory = dirname(fileURLToPath(import.meta.url));
const outputDirectory = resolve(docsDirectory, "dist");
const katexDistribution = resolve(docsDirectory, "node_modules/katex/dist");
const katexOutput = resolve(outputDirectory, "katex");

const protectedBlocks = /<(pre|code|script|style|textarea)\b[^>]*>[\s\S]*?<\/\1>/gi;
const doxygenFormula = /FUGACITYKATEX(DISPLAY|INLINE)([A-Za-z0-9_-]+)END/g;
const displayFormula = /\$\$([\s\S]*?)\$\$/g;
const inlineFormula = /(?<!\\)\$(?!\$)([\s\S]*?)(?<!\\)\$(?!\$)/g;

function decodeEntities(value) {
  return value
    .replace(/&(?:#x([0-9a-f]+)|#([0-9]+));/gi,
      (_entry, hex, decimal) => String.fromCodePoint(
        Number.parseInt(hex ?? decimal, hex ? 16 : 10),
      ))
    .replaceAll("&lt;", "<")
    .replaceAll("&gt;", ">")
    .replaceAll("&quot;", '"')
    .replaceAll("&#39;", "'")
    .replaceAll("&amp;", "&");
}

function markdownHtmlToTex(fragment) {
  return decodeEntities(fragment
    .replace(/<em>/g, "_")
    .replace(/<\/em>/g, "_")
    .replace(/<br\s*\/?>/gi, "\n")
    .replace(/<[^>]+>/g, ""))
    .trim();
}

function render(tex, displayMode, source) {
  try {
    return katex.renderToString(tex, {
      displayMode,
      output: "htmlAndMathml",
      strict: "error",
      throwOnError: true,
      trust: false,
    });
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    throw new Error(`Invalid TeX in ${source}: ${message}\n${tex}`);
  }
}

async function htmlFiles(directory) {
  const files = [];
  for (const entry of await readdir(directory, { withFileTypes: true })) {
    const path = resolve(directory, entry.name);
    if (entry.isDirectory()) {
      files.push(...await htmlFiles(path));
    } else if (entry.isFile() && entry.name.endsWith(".html")) {
      files.push(path);
    }
  }
  return files;
}

await mkdir(katexOutput, { recursive: true });
await cp(resolve(katexDistribution, "katex.min.css"),
  resolve(katexOutput, "katex.min.css"));
await cp(resolve(katexDistribution, "fonts"), resolve(katexOutput, "fonts"), {
  recursive: true,
});

let displayCount = 0;
let inlineCount = 0;
const files = await htmlFiles(outputDirectory);

for (const path of files) {
  const protectedContent = [];
  let html = await readFile(path, "utf8");
  html = html.replace(protectedBlocks, (block) => {
    const index = protectedContent.push(block) - 1;
    return `FUGACITYPROTECTEDBLOCK${index}END`;
  });

  html = html.replace(doxygenFormula, (_entry, mode, payload) => {
    const display = mode === "DISPLAY";
    if (display) displayCount += 1;
    else inlineCount += 1;
    const tex = decodeEntities(Buffer.from(payload, "base64url").toString("utf8"));
    return render(tex, display, path);
  });
  html = html.replace(displayFormula, (_entry, fragment) => {
    displayCount += 1;
    return render(markdownHtmlToTex(fragment), true, path);
  });
  html = html.replace(inlineFormula, (_entry, fragment) => {
    inlineCount += 1;
    return render(markdownHtmlToTex(fragment), false, path);
  });

  if (doxygenFormula.test(html) || displayFormula.test(html) || inlineFormula.test(html)) {
    throw new Error(`Unrendered TeX delimiter or marker remains in ${path}`);
  }

  const relativeRoot = relative(dirname(path), outputDirectory).replaceAll("\\", "/");
  const katexHref = `${relativeRoot ? `${relativeRoot}/` : ""}katex/katex.min.css`;
  html = html.replace("</head>",
    `<link rel="stylesheet" href="${katexHref}"/></head>`);
  html = html.replace(/FUGACITYPROTECTEDBLOCK([0-9]+)END/g,
    (_entry, index) => protectedContent[Number(index)]);
  await writeFile(path, html);
}

process.stdout.write(
  `Rendered ${displayCount} display and ${inlineCount} inline formulas with KaTeX `
    + `across ${files.length} HTML pages.\n`,
);
