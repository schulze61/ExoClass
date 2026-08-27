"""
nb_toc.py — Jupyter notebook TOC utilities.

Usage
-----
from nb_toc import toc, back_to_toc, setup

# REQUIRED: Run this once in its own cell before toc() and back_to_toc():
setup()

# In your TOC cell (code cell, placed just below ## Table of Contents):
toc([
    "Introduction",
    "Data and Methods",
    "Results",
    "Discussion",
    "Conclusion",
])

# At the bottom of each section (code cell):
back_to_toc()
"""

from IPython.display import Javascript, HTML, display
import uuid


def setup():
    """
    Run this once in its own cell before calling toc() or back_to_toc().
    Registers the core scrolling logic used by all TOC buttons.
    """
    display(Javascript("""
(function() {
  if (window._jptocReady) { return; }
  window._jptocReady = true;

  window._jptocScrollTo = function(headingText, tocId) {
    // Find the notebook panel containing the TOC with this ID
    var tocEl = document.getElementById(tocId);
    var panel = tocEl ? tocEl.closest('.jp-NotebookPanel') : null;
    var scope = panel || document;

    function findHeading() {
      var headings = scope.querySelectorAll('h1,h2,h3,h4,h5,h6');
      for (var i = 0; i < headings.length; i++) {
        if (headings[i].textContent.replace(/\u00b6/g, '').trim() === headingText) {
          return headings[i];
        }
      }
      return null;
    }

    var el = findHeading();
    if (el) {
      el.scrollIntoView();
      return;
    }

    var outer = scope.querySelector('.jp-WindowedPanel-outer');
    var inner = scope.querySelector('.jp-WindowedPanel-inner');
    if (!outer || !inner) { return; }

    var step           = outer.clientHeight * 0.8;
    var direction      = 1;
    var passesLeft     = 3;
    var pos            = outer.scrollTop;
    var savedScrollTop = outer.scrollTop;

    outer.style.visibility = 'hidden';

    function restore(targetEl) {
      outer.style.visibility = '';
      if (targetEl) {
        targetEl.scrollIntoView();
      } else {
        outer.scrollTop = savedScrollTop;
      }
    }

    function sweep() {
      var el = findHeading();
      if (el) {
        restore(el);
        return;
      }

      var totalHeight = inner.offsetHeight;
      var viewHeight  = outer.clientHeight;
      var maxScroll   = totalHeight - viewHeight;

      pos += step * direction;

      if (pos >= maxScroll) {
        pos = maxScroll;
        direction = -1;
        passesLeft--;
      }

      if (pos <= 0) {
        pos = 0;
        direction = 1;
        passesLeft--;
      }

      if (passesLeft <= 0) {
        restore(null);
        return;
      }

      outer.scrollTop = pos;
      setTimeout(sweep, 150);
    }

    setTimeout(sweep, 150);
  };

  document.addEventListener('click', function(e) {
    var btn = e.target.closest('.jptoc-btn, .jptoc-back');
    if (!btn) { return; }
    window._jptocScrollTo(
      btn.getAttribute('data-heading'),
      btn.getAttribute('data-toc-id')
    );
  });
})();
"""))


def toc(
    sections: list[str],
    title: str = "Table of Contents",
) -> HTML:
    """
    Render a clickable Table of Contents.

    Parameters
    ----------
    sections : list of str
        Heading text for each section, exactly as written in the ## cells.
    title : str
        Display title shown above the TOC list.

    Returns
    -------
    IPython.display.HTML
        Call as the last expression in a code cell to display.
    """
    toc_id = "jptoc-" + uuid.uuid4().hex[:8]

    items_html = "\n".join(
        f'    <li>'
        f'<button class="jptoc-btn" data-heading="{s}" data-toc-id="{toc_id}" '
        f'style="background:none;border:none;padding:0;color:#0066cc;'
        f'cursor:pointer;font-size:1em;text-align:left;">'
        f'{s}'
        f'</button></li>'
        for s in sections
    )

    return HTML(f"""\
<div id="{toc_id}" style="background:#f8f8f8;border:1px solid #ddd;border-radius:6px;
            padding:14px 18px;display:inline-block;min-width:220px">
  <div style="font-weight:bold;font-size:1.05em;margin-bottom:8px">{title}</div>
  <ol style="margin:0;padding-left:18px;line-height:1.9">
{items_html}
  </ol>
</div>""")


def back_to_toc(toc_heading: str = "Table of Contents") -> HTML:
    """
    Render an upward arrow button that scrolls back to the TOC heading.

    Parameters
    ----------
    toc_heading : str
        The heading text of the TOC cell. Must match the ## heading exactly.
        Defaults to "Table of Contents".

    Returns
    -------
    IPython.display.HTML
        Call as the last expression in a code cell to display.
    """
    # Use a placeholder toc_id since back_to_toc just needs to find a heading
    # The panel scoping will be determined by the closest NotebookPanel
    return HTML(
        f'<button class="jptoc-back" data-heading="{toc_heading}" data-toc-id="jptoc-back" '
        f'style="background:none;border:none;padding:0;color:#555;'
        f'cursor:pointer;font-size:0.9em;">&#8593; Back to TOC</button>'
    )
