import assert from "node:assert/strict";
import { mkdtemp, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { resolve } from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const fixture = await mkdtemp(resolve(tmpdir(), "fugacity-doxygen-filter-"));
try {
  await writeFile(resolve(fixture, "index.xml"), "<doxygenindex></doxygenindex>\n");
  await writeFile(resolve(fixture, "orphan.xml"),
    '<doxygen><para><ref kindref="compound" refid="missing">missing</ref></para></doxygen>\n');
  const filter = resolve(fileURLToPath(new URL(".", import.meta.url)), "filter-doxygen-xml.mjs");
  const result = spawnSync(process.execPath, [filter, fixture], { encoding: "utf8" });
  assert.notEqual(result.status, 0);
  assert.match(result.stderr, /1 orphan references; refusing to build/);
} finally {
  await rm(fixture, { recursive: true, force: true });
}

process.stdout.write("Doxygen XML orphan-reference test passed.\n");
