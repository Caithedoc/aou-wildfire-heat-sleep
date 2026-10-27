#!/bin/bash
# PRISM 4km 日最低温 2016-2023 批量下载 (与 download_prism.sh 同构, tmean->tmin)
# https://ftp.prism.oregonstate.edu/time_series/us/an/4km/tmin/daily/YYYY/prism_tmin_us_25m_YYYYMMDD.zip
set -u
DEST="$(dirname "$0")/raw/prism_tmin"
mkdir -p "$DEST"
cd "$DEST"
: > fail.log
python3 - <<'PY' > /tmp/prism_tmin_dates.txt
from datetime import date, timedelta
d, end = date(2016, 1, 1), date(2023, 12, 31)
while d <= end:
    print(d.strftime("%Y%m%d"))
    d += timedelta(days=1)
PY
N=0
while read -r ymd; do
  y=${ymd:0:4}
  f="prism_tmin_us_25m_${ymd}.zip"
  if [ -s "$f" ]; then continue; fi
  url="https://ftp.prism.oregonstate.edu/time_series/us/an/4km/tmin/daily/${y}/${f}"
  curl -sf --retry 2 --retry-delay 3 -o "$f" "$url" || echo "$ymd" >> fail.log
  N=$((N+1))
  if [ $((N % 200)) -eq 0 ]; then echo "downloaded $N files"; fi
  sleep 0.15
done < /tmp/prism_tmin_dates.txt
# 补漏一轮
if [ -s fail.log ]; then
  cp fail.log fail_round1.log
  : > fail.log
  while read -r ymd; do
    y=${ymd:0:4}
    f="prism_tmin_us_25m_${ymd}.zip"
    url="https://ftp.prism.oregonstate.edu/time_series/us/an/4km/tmin/daily/${y}/${f}"
    curl -sf --retry 3 --retry-delay 5 -o "$f" "$url" || echo "$ymd" >> fail.log
    sleep 0.3
  done < fail_round1.log
fi
echo "DONE. files: $(ls prism_tmin_us_25m_*.zip 2>/dev/null | wc -l), failures: $(wc -l < fail.log)"
