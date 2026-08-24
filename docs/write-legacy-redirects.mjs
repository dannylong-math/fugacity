import { access, mkdir, writeFile } from "node:fs/promises";
import { dirname, relative, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const docsDirectory = dirname(fileURLToPath(import.meta.url));
const outputDirectory = resolve(docsDirectory, "dist");
const redirects = new Map([
  // Sourcey 3.6.5's Doxygen navigation targets api/index.html even when its
  // flat-URL index artifact is api.html.
  ["api/index.html", "api.html"],
  ["getting_started.html", "getting-started.html"],
  ["implementing_a_new_eos.html", "implementing-a-new-eos.html"],
  ["api/concepts.html", "concepts.html"],
]);

// tutorial.html is both the legacy path and Sourcey's canonical flat-URL path.
await access(resolve(outputDirectory, "tutorial.html"));

for (const [legacyPath, canonicalPath] of redirects) {
  const destination = resolve(outputDirectory, legacyPath);
  const canonical = resolve(outputDirectory, canonicalPath);
  await access(canonical);
  await mkdir(dirname(destination), { recursive: true });
  const href = relative(dirname(destination), canonical).replaceAll("\\", "/");
  const html = `<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta http-equiv="refresh" content="0; url=${href}">
<link rel="canonical" href="${href}"><title>Documentation moved</title></head>
<body><p>This page moved to <a href="${href}">${href}</a>.</p></body></html>
`;
  await writeFile(destination, html);
}

process.stdout.write(`Wrote ${redirects.size} documentation compatibility redirects.\n`);
