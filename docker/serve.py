#!/usr/bin/env python3
"""Minimal, zero-dependency static server for the admissible project site.

In keeping with the tool's zero-runtime-dependency stance, the website is
served by the standard library alone — no framework, no wheels. It hosts the
bundled ``web/`` directory and answers a cheap health check at ``/up`` (and
``/healthz``), which is what the platform's deploy gate and monitor poll.
"""
from __future__ import annotations

import os
from functools import partial
from http.server import HTTPServer, SimpleHTTPRequestHandler

HERE = os.path.dirname(os.path.abspath(__file__))
WEB_ROOT = os.path.abspath(os.path.join(HERE, os.pardir, "web"))
PORT = int(os.environ.get("PORT", "8000"))


class Handler(SimpleHTTPRequestHandler):
    """Serve static files, with a plain-text health endpoint bolted on."""

    def _health(self) -> None:
        body = b"ok\n"
        self.send_response(200)
        self.send_header("Content-Type", "text/plain; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def do_GET(self) -> None:  # noqa: N802 - name fixed by BaseHTTPRequestHandler
        if self.path in ("/up", "/healthz"):
            self._health()
            return
        super().do_GET()

    def do_HEAD(self) -> None:  # noqa: N802 - name fixed by BaseHTTPRequestHandler
        if self.path in ("/up", "/healthz"):
            self._health()
            return
        super().do_HEAD()

    def log_message(self, fmt: str, *args) -> None:
        # Skip the per-request noise; the platform logs traffic at the edge.
        pass


def main() -> None:
    handler = partial(Handler, directory=WEB_ROOT)
    httpd = HTTPServer(("0.0.0.0", PORT), handler)
    print(f"admissible site: serving {WEB_ROOT} on 0.0.0.0:{PORT}", flush=True)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        httpd.shutdown()


if __name__ == "__main__":
    main()
