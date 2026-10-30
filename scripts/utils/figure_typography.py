"""Typography for figures exported at their final manuscript dimensions."""
from matplotlib.text import Text


def readable_text(fig, minimum=8.0):
    """Keep labels readable; call before saving, at the intended physical size."""
    # Materialise tick labels before setting their sizes. Constrained layouts
    # then recompute spacing during export using the enlarged text extents.
    fig.canvas.draw()
    for text in fig.findobj(Text):
        if text.get_text().strip():
            text.set_fontsize(max(minimum, text.get_fontsize()))
