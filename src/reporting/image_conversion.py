import base64
import re
from pathlib import Path

def inline_base64_images(html: str, base_dir=".") -> str:
    pattern = re.compile(r'<img\s+src="([^"]+)"(?:\s+width="(\d+)")?(?:\s+height="(\d+)")?\s*>')

    def repl(match):
        src, width, height = match.groups()
        image_path = Path(base_dir) / src
        if not image_path.exists():
            print(f"Warning: Could not find image {image_path}")
            return match.group(0)

        ext = image_path.suffix.lower().lstrip('.')
        mime = "image/jpeg" if ext in ['jpg', 'jpeg'] else f"image/{ext}"

        with open(image_path, "rb") as f:
            encoded = base64.b64encode(f.read()).decode("utf-8")

        style_parts = []
        if width:
            style_parts.append(f"width:{width}px")
        if height:
            style_parts.append(f"height:{height}px")
        style_attr = f' style="{"; ".join(style_parts)}"' if style_parts else ""

        return f'<img src="data:{mime};base64,{encoded}"{style_attr}>'

    return pattern.sub(repl, html)
