import { cp, mkdir, readdir, readFile, writeFile } from "node:fs/promises";
import { dirname, relative, resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

import katex from "katex";

const docsDirectory = dirname(fileURLToPath(import.meta.url));
const outputDirectory = resolve(docsDirectory, "dist");
const katexDistribution = resolve(docsDirectory, "node_modules/katex/dist");
const katexOutput = resolve(outputDirectory, "katex");

const protectedBlocks = /<(pre|code|script|style|textarea)\b[^>]*>[\s\S]*?<\/\1>/gi;
const doxygenFormula = /FUGACITYKATEX:(DISPLAY|INLINE):([A-Za-z0-9_-]+):FUGACITYEND/g;
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

function restoreTags(fragment, tags) {
  return fragment.replace(/FUGACITYPROTECTEDTAG([0-9]+)END/g,
    (_entry, index) => tags[Number(index)]);
}

function markdownHtmlToTex(fragment, tags) {
  return decodeEntities(restoreTags(fragment, tags)
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

// Quote-aware tag protection lets formulas span inline markup while ensuring
// dollar signs in HTML attributes are never interpreted as TeX delimiters.
function protectTags(html) {
  const tags = [];
  let protectedHtml = "";
  let cursor = 0;
  while (cursor < html.length) {
    if (html[cursor] !== "<") {
      protectedHtml += html[cursor];
      cursor += 1;
      continue;
    }
    let end = cursor + 1;
    let quote = null;
    for (; end < html.length; end += 1) {
      const character = html[end];
      if (quote) {
        if (character === quote) quote = null;
      } else if (character === '"' || character === "'") {
        quote = character;
      } else if (character === ">") {
        break;
      }
    }
    if (end === html.length) throw new Error("Malformed HTML: unterminated tag");
    const tag = html.slice(cursor, end + 1);
    const index = tags.push(tag) - 1;
    protectedHtml += `FUGACITYPROTECTEDTAG${index}END`;
    cursor = end + 1;
  }
  return { protectedHtml, tags };
}

function firstUnescapedDollar(value) {
  for (let index = 0; index < value.length; index += 1) {
    if (value[index] !== "$") continue;
    let backslashes = 0;
    for (let before = index - 1; before >= 0 && value[before] === "\\"; before -= 1) {
      backslashes += 1;
    }
    if (backslashes % 2 === 0) return index;
  }
  return -1;
}

export function renderHtmlMath(input, source = "HTML input") {
  const protectedContent = [];
  let html = input.replace(protectedBlocks, (block) => {
    const index = protectedContent.push(block) - 1;
    return `FUGACITYPROTECTEDBLOCK${index}END`;
  });
  const { protectedHtml, tags } = protectTags(html);
  html = protectedHtml;

  let displayCount = 0;
  let inlineCount = 0;
  html = html.replace(doxygenFormula, (_entry, mode, payload) => {
    const display = mode === "DISPLAY";
    if (display) displayCount += 1;
    else inlineCount += 1;
    const tex = decodeEntities(Buffer.from(payload, "base64url").toString("utf8"));
    return render(tex, display, source);
  });
  html = html.replace(displayFormula, (_entry, fragment) => {
    displayCount += 1;
    return render(markdownHtmlToTex(fragment, tags), true, source);
  });
  html = html.replace(inlineFormula, (_entry, fragment) => {
    inlineCount += 1;
    return render(markdownHtmlToTex(fragment, tags), false, source);
  });

  if (html.includes("FUGACITYKATEX:")) {
    throw new Error(`Unrendered Doxygen formula marker remains in ${source}`);
  }
  const unmatched = firstUnescapedDollar(html);
  if (unmatched >= 0) {
    const context = html.slice(Math.max(0, unmatched - 30), unmatched + 31);
    throw new Error(`Unmatched TeX delimiter remains in ${source}: ${context}`);
  }

  html = restoreTags(html, tags);
  html = html.replace(/FUGACITYPROTECTEDBLOCK([0-9]+)END/g,
    (_entry, index) => protectedContent[Number(index)]);
  return { html, displayCount, inlineCount };
}

async function htmlFiles(directory) {
  const files = [];
  for (const entry of await readdir(directory, { withFileTypes: true })) {
    const path = resolve(directory, entry.name);
    if (entry.isDirectory()) files.push(...await htmlFiles(path));
    else if (entry.isFile() && entry.name.endsWith(".html")) files.push(path);
  }
  return files;
}

async function main() {
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
    const input = await readFile(path, "utf8");
    const rendered = renderHtmlMath(input, path);
    displayCount += rendered.displayCount;
    inlineCount += rendered.inlineCount;
    const relativeRoot = relative(dirname(path), outputDirectory).replaceAll("\\", "/");
    const katexHref = `${relativeRoot ? `${relativeRoot}/` : ""}katex/katex.min.css`;
    const html = rendered.html.replace("</head>",
      `<link rel="stylesheet" href="${katexHref}"/></head>`);
    await writeFile(path, html);
  }

  process.stdout.write(
    `Rendered ${displayCount} display and ${inlineCount} inline formulas with KaTeX `
      + `across ${files.length} HTML pages.\n`,
  );
}

const invokedDirectly = process.argv[1]
  && import.meta.url === pathToFileURL(resolve(process.argv[1])).href;
if (invokedDirectly) await main();
