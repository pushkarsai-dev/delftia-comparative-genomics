"""Regression tests binding saved exports to the figure's actual artists."""

import contextlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from PIL import Image, ImageDraw

PATH = Path(__file__).resolve().parents[1] / "scripts/validate_figure1.py"
SPEC = importlib.util.spec_from_file_location("figure_validator", PATH)
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


class ExportIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="figure1-export-test-")
        self.directory = Path(self.temporary.name)
        self.pdf = self.directory / "figure.pdf"
        self.png = self.directory / "figure.png"
        self.figure = Figure(figsize=(1.5, 1.5))
        self.axis = self.figure.subplots()
        self.axis.plot([0, 1], [0, 1], color="blue", linewidth=4)
        self.figure.savefig(self.pdf, format="pdf")
        self.figure.savefig(self.png, format="png", dpi=600)

    def tearDown(self):
        self.figure.clear()
        self.temporary.cleanup()
        VALIDATOR.FAIL.clear()

    def validate(self):
        VALIDATOR.FAIL.clear()
        with contextlib.redirect_stdout(io.StringIO()):
            VALIDATOR.export_correspondence({"fig": self.figure}, self.pdf, self.png)
        return list(VALIDATOR.FAIL)

    def test_unchanged_exports_pass(self):
        self.assertEqual(self.validate(), [])

    def test_single_changed_png_pixel_is_rejected(self):
        with Image.open(self.png) as source:
            image = source.convert("RGBA")
        image.putpixel((0, 0), (0, 0, 0, 255))
        image.save(self.png)
        failures = self.validate()
        self.assertEqual(failures, ["all saved PNG pixels equal the current figure export"])

    def test_blank_png_region_is_rejected(self):
        with Image.open(self.png) as source:
            image = source.convert("RGBA")
        ImageDraw.Draw(image).rectangle((0, 0, image.width - 1, int(image.height * .79)), fill="white")
        image.save(self.png)
        self.assertIn("all saved PNG pixels equal the current figure export", self.validate())

    def test_changed_saved_pdf_is_rejected(self):
        altered = Figure(figsize=(1.5, 1.5))
        altered.subplots().plot([0, 1], [1, 0], color="red", linewidth=4)
        altered.savefig(self.pdf, format="pdf")
        altered.clear()
        failures = self.validate()
        self.assertEqual(failures, ["all rendered PDF pixels equal current export, page 1 (200 dpi)"])


if __name__ == "__main__":
    unittest.main()
