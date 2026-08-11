#!/usr/bin/env python3
"""Serve field/console.html on localhost (rebuilds first).

  python scripts/serve_field_console.py
  python scripts/serve_field_console.py --port 8765
"""
from __future__ import annotations

import argparse
import http.server
import socketserver
import sys
import webbrowser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_field_console import main as build_main  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--no-open", action="store_true")
    args = ap.parse_args()

    build_main([])
    field = ROOT / "field"
    if not (field / "console.html").exists():
        print("console.html missing after build", file=sys.stderr)
        return 1

    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *a, **k):
            super().__init__(*a, directory=str(field), **k)

        def log_message(self, fmt, *args):
            sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    url = f"http://127.0.0.1:{args.port}/console.html"
    print(f"Serving {field}")
    print(f"Open {url}")
    print("Ctrl+C to stop")
    if not args.no_open:
        webbrowser.open(url)
    with socketserver.TCPServer(("127.0.0.1", args.port), Handler) as httpd:
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nstopped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
