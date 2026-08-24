import assert from "node:assert/strict";

import { renderHtmlMath } from "./render-math.mjs";

const inline = renderHtmlMath("<p>$x^2$</p>", "inline fixture");
assert.match(inline.html, /class="katex"/);
assert.equal(inline.inlineCount, 1);

const display = renderHtmlMath("<p>$$x^2$$</p>", "display fixture");
assert.match(display.html, /class="katex-display"/);
assert.equal(display.displayCount, 1);

const escapedDollar = renderHtmlMath("<p>$\\$5$</p>", "escaped dollar fixture");
assert.match(escapedDollar.html, /class="katex"/);
assert.match(escapedDollar.html, /\$5/);

const textDollar = renderHtmlMath(
  "<p>$\\text{price \\$5}$</p>",
  "text dollar fixture",
);
assert.match(textDollar.html, /class="katex"/);
assert.match(textDollar.html, /price(?: |\u00a0)\$5/);

assert.throws(
  () => renderHtmlMath("<p>$x^2</p>", "unmatched inline fixture"),
  /Unmatched TeX delimiter/,
);
assert.throws(
  () => renderHtmlMath("<p>$$x^2</p>", "unmatched display fixture"),
  /Unmatched TeX delimiter/,
);

const protectedInput = '<a title="$attribute">cost</a><pre>$x^2</pre><code>$$</code>';
assert.equal(renderHtmlMath(protectedInput, "protected fixture").html, protectedInput);

const formerSentinels = [
  "<p>FUGACITYPROTECTEDTAG0END FUGACITYPROTECTEDBLOCK0END FUGACITYTEMP</p>",
  "<pre>protected</pre>",
].join("");
assert.equal(renderHtmlMath(formerSentinels, "sentinel fixture").html, formerSentinels);

assert.throws(
  () => renderHtmlMath("<p>$\\notacommand$</p>", "invalid TeX fixture"),
  /Invalid TeX/,
);

process.stdout.write("Math renderer adversarial tests passed.\n");
