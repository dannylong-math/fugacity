import assert from "node:assert/strict";

import { renderHtmlMath } from "./render-math.mjs";

const inline = renderHtmlMath("<p>$x^2$</p>", "inline fixture");
assert.match(inline.html, /class="katex"/);
assert.equal(inline.inlineCount, 1);

const display = renderHtmlMath("<p>$$x^2$$</p>", "display fixture");
assert.match(display.html, /class="katex-display"/);
assert.equal(display.displayCount, 1);

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

assert.throws(
  () => renderHtmlMath("<p>$\\notacommand$</p>", "invalid TeX fixture"),
  /Invalid TeX/,
);

process.stdout.write("Math renderer adversarial tests passed.\n");
