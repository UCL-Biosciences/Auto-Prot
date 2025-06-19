import base64
import re
from pathlib import Path
from PIL import Image

def inline_base64_images(html: str, base_dir=".") -> str:
    pattern = re.compile(r'<img\s+src="([^"]+)"(?:\s+width="(\d+)")?(?:\s+height="(\d+)")?\s*>')

    def repl(match):
        src, width_attr, height_attr = match.groups()
        image_path = Path(base_dir) / src
        if not image_path.exists():
            print(f"Warning: Could not find image {image_path}")
            return match.group(0)

        ext = image_path.suffix.lower().lstrip('.')
        mime = "image/jpeg" if ext in ['jpg', 'jpeg'] else f"image/{ext}"

        with open(image_path, "rb") as f:
            encoded = base64.b64encode(f.read()).decode("utf-8")

        # Check actual image size
        try:
            with Image.open(image_path) as img:
                actual_width, actual_height = img.size
        except Exception as e:
            print(f"Error opening image for size check: {image_path} — {e}")
            actual_width = None

        # Decide styling
        style_parts = []
        if actual_width and actual_width > 900:
            style_parts.append("width:900px")  # Let browser scale height accordingly
        style_attr = f' style="{"; ".join(style_parts)}"' if style_parts else ""

        return f'<img src="data:{mime};base64,{encoded}"{style_attr}>'

    return pattern.sub(repl, html)

