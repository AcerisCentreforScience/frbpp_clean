"""Pull retained GCN Kafka history into data/raw/gcn.

Usage: python scripts/fetch_gcn_history.py
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import sys
import time
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from gcn_live import make_oauth_cb  # noqa: E402

CFG = json.loads((ROOT / "config" / "gcn_kafka.local.json").read_text(encoding="utf-8"))
TOPICS = ["gcn.notices.chime.frb", "gcn.notices.dsa110.frb"]


def main() -> int:
    from gcn_kafka import Consumer

    consumer = Consumer(config={"oauth_cb": make_oauth_cb(CFG),
                                "auto.offset.reset": "earliest",
                                "group.id": "hist-" + uuid.uuid4().hex[:8]},
                        client_id=CFG["client_id"], client_secret=CFG["client_secret"])
    consumer.subscribe(TOPICS)
    print("订阅:", TOPICS)

    n, t0, last = 0, time.time(), time.time()
    ts_seen: list[int] = []
    lines: list[str] = []
    while time.time() - t0 < 300:
        msgs = consumer.consume(timeout=1.0, num_messages=500)
        if not msgs:
            if time.time() - last > 20:
                print(f"静默 20 s（已收到 {n} 条），视为拉尽")
                break
            continue
        last = time.time()
        for m in msgs:
            if m.error():
                continue
            n += 1
            ts = m.timestamp()[1] if m.timestamp() else None
            if ts:
                ts_seen.append(ts)
            lines.append(json.dumps(
                {"topic": m.topic(), "partition": m.partition(), "offset": m.offset(),
                 "kafka_ts_ms": ts, "value": m.value().decode("utf-8", "replace")[:4000]},
                ensure_ascii=False))
            if n <= 15:
                print(f"[{m.topic()}] off={m.offset()} {lines[-1][:200]}")
    consumer.close()

    if not ts_seen:
        print("未收到任何消息；不落盘")
        return 1
    lo = dt.datetime.fromtimestamp(min(ts_seen) / 1000, dt.timezone.utc).strftime("%Y%m%d")
    hi = dt.datetime.fromtimestamp(max(ts_seen) / 1000, dt.timezone.utc).strftime("%Y%m%d")
    outdir = ROOT / "data" / "raw" / "gcn"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"history_{lo}_{hi}.jsonl"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    sha = hashlib.sha256(out.read_bytes()).hexdigest()
    print(f"总消息 {n} → {out.relative_to(ROOT)}  sha256={sha}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())