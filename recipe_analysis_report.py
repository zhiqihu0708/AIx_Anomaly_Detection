"""
recipe_analysis_report.py
==========================
Unified recipe trace + step importance analysis.  Reads a sensor time-series
CSV and produces a single presentation-style HTML report covering:

  1. Features & Logic Overview
     Explains every algorithm and design decision in the tool.

  2. All-Runs Status Overview
     Every run shown in one table (Normal / Anomalous / Aborted).
     Baseline envelope chart for the primary fault sensor overlaid with all
     runs.  One representative z-score chart for the worst anomalous run.

  3. Lead / Lag Cascade Detail
     Sensor divergence timeline + alarm interleave for each anomalous run.

  4. Step Importance Ranking
     Recipe steps ranked by cross-run coefficient of variation  --  identifies
     which step carries the most diagnostic signal.

  5. Per-Step Trend Charts, Highfliers & Correlation Heatmaps
     Top-N steps: per-run mean trend, highflier classification, Pearson
     correlation matrix (aborted + data-artifact runs excluded).

  6. Baseline Envelope Traces
     Median ± 3×MAD green band vs individual run traces for all deviating
     sensors.

  7. Raw Sensor Traces & SP Deviation
     Full run overlays for cascade sensors, actual vs setpoint charts.

  8. Sensor Deviation Classification
     Genuine Fault vs Change-in-SP table.

Zero external dependencies  --  pure Python 3.6+.

Usage:
    python recipe_analysis_report.py
    → prompts for CSV file, writes <stem>_combined_report.html
"""

import csv, bisect, math, os, re, sys, msvcrt
from collections import defaultdict

# Ensure Unicode output works on Windows console
if hasattr(sys.stdout, "reconfigure"):
    try: sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception: pass

# ══════════════════════════════════════════════════════════════════════════════
# TUNABLES
# ══════════════════════════════════════════════════════════════════════════════

# ── recipe-trace tunables ─────────────────────────────────────────────────────
SP_OFF_FRAC_THRESHOLD  = 0.05
SP_TRACKING_THRESHOLD  = 0.05
KNOWN_ANOMALOUS: set   = set()
SENSOR_DIVERGE_SIGMA   = 3.0
OUTLIER_RUN_ZSCORE     = 2.5
MIN_CONSEC_DIVERGE     = 3
BASELINE_Z_THRESH      = 1.5
CHART_HEIGHT           = 480
MIN_PX_PER_STEP        = 50
MAX_CHART_WIDTH        = 8000
ZSCORE_WIDTH           = 1400
ZSCORE_HEIGHT          = 520

# ── step-importance tunables ──────────────────────────────────────────────────
SKIP_STEPS             = {"-1", "0"}   # bookkeeping steps never scored
MIN_STEP_ROWS          = 500
MIN_RUNS_FOR_SCORE     = 5
TOP_N_STEPS            = 7
MIN_CV_FOR_CHART       = 0.005
MIN_STD_FRAC_OF_RANGE  = 0.01
MIN_RNG_FRAC_OF_MEAN   = 0.03  # run-means range must be >= 3% of step mean to show trend chart
MIN_ACTIVE_FRAC        = 0.05
ABORT_ROW_FRAC         = 0.75
FLOW_MIN_RANGE_SCCM    = 10.0
MFC_RAMP_RATE_SCCM_S   = 5.0   # sccm/s -- if median flow changes faster than this, it's a ramp
MFC_RAMP_WINDOW_S      = 10.0  # seconds -- ramp window padding on each side
MFC_OVERSHOOT_FRAC     = 0.20  # flag ramp rows where flow > sp_target * (1 + this fraction)
SCORE_OUTLIER_RATIO    = 10.0  # if one sensor's max score is > this x median of others, segregate
HIGHFLIER_IQR_MULT     = 3.0
STEP_NOISE_RATIO       = 0.20  # within-run std / step range > this => noisy in that step
STEP_INACTIVE_FRAC     = 0.015 # abs(median run-mean) / global max < this => inactive in that step
STEP_NOISY_STEP_FRAC   = 0.80  # if sensor noisy/inactive in > this fraction of steps, exclude globally
N_SENSORS_FOR_REAL     = 2
STEP_ROWS_MIN_FRAC     = 0.80
TREND_CHART_W          = 940
TREND_CHART_H          = 280
TREND_MARGIN           = dict(top=44, right=24, bottom=62, left=74)

# ── shared colors ─────────────────────────────────────────────────────────────
RUN_COLORS = [
    "#1565C0","#2E7D32","#6A1B9A","#00838F","#4E342E",
    "#37474F","#558B2F","#0277BD","#AD1457","#E65100",
]
LEAD_COLOR       = "#D32F2F"
ABORT_COLOR      = "#E65100"
HIGHFLIER_REAL   = "#B71C1C"
HIGHFLIER_ART    = "#FBC02D"
TOOL_COLORS      = RUN_COLORS

def _cascade_color(ci, n):
    if n <= 1: return "#E65100"
    hue = int(40 * (1 - ci / max(n - 1, 1)))
    return f"hsl({hue},90%,40%)"

# ── keyword lists (schema discovery) ─────────────────────────────────────────
NUMERIC_FRAC_THRESHOLD = 0.70

TIME_KEYWORDS = [
    "time", "timestamp", "ts", "elapsed", "seconds", "sec", "t_s",
    "time_s", "datetime", "date", "epoch",
]
SKIP_KEYWORDS = ["comment", "note", "description", "remark", "memo"]
ID_KEYWORDS   = [
    "tool", "chamber", "slot", "station", "port", "module", "unit",
    "position", "location", "zone", "head", "chuck", "carrier",
    "cassette", "boat", "rack", "nest", "pocket", "lane", "arm",
    "id", "index", "idx", "num", "number", "no", "seq",
    "run", "lot", "wafer", "batch", "recipe", "program",
    "site", "die", "device", "equipment", "eqp",
]
LABEL_KEYWORDS = [
    "alarm", "flag", "alert", "status", "code", "phase", "mode",
    "state", "recipe", "wafer", "lot", "step", "scenario", "type",
    "category", "class", "label", "result", "pass", "fail",
    "severity", "priority", "level",
]
# Keywords that suggest a column carries alarm/event text
ALARM_KEYWORDS  = ["alarm", "alarmcode", "alarmstatus", "alert", "fault"]
EVENT_KEYWORDS  = ["event", "eventname", "eventtype", "eventsource",
                   "eventdescription", "eventid"]
# Keywords for step columns
STEP_NAME_KW    = ["step_name", "stepname", "step_description",
                   "recipe_step", "process_step", "step_label"]
STEP_NUM_KW     = ["step_number", "stepnumber", "step_num", "step_no",
                   "step_id", "step_index", "step_seq"]
STEP_ELAPSED_KW = ["recipe_step_elapsed", "step_elapsed", "step_time",
                   "recipe_step_time", "step_duration",
                   "recipe_step_elapsed_time"]
STEP_SEQ_KW     = ["calcstepseq", "stepseq", "calc_step_seq",
                   "step_sequence", "step_seq"]

# Keywords that identify event-counter / tally columns that should be
# EXCLUDED from anomaly scoring and lead/lag detection.  These columns are
# discrete integer tallies (arc event counts, alarm counts, fault counts),
# not continuous sensor signals; including them skews scoring and produces
# spurious "lead sensor" assignments.
EVT_COUNT_KW = [
    "evt_count", "evtcount", "event_count", "eventcount",
    "arc_evt",   "arcevt",   "arc_count",   "arccount",
    "fault_count","faultcount",
    "alarm_count","alarmcount",
    "error_count","errorcount",
    "tally",     "cumcount",  "cumulative_count",
    "ioname_arc",            # e.g. HighFreqRF_IONAME_ARC_EVT_COUNT
]

# User-supplied exclusions -- populated at runtime via prompt in main().
EXCLUDED_SENSORS: set = set()


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1  --  SHARED HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def _short(s, maxlen=30):
    """Shorten a sensor name for chart/table labels."""
    if len(s) <= maxlen:
        return s
    for sep in ('/', '\\', '.', '_'):
        idx = s.rfind(sep)
        if idx != -1 and len(s) - idx - 1 >= 4:
            return s[idx + 1:]
    return s[-maxlen:]

def _iqr(vals):
    """Interquartile range (Q3 - Q1) of a numeric sequence."""
    s = sorted(vals)
    n = len(s)
    mid = n // 2
    lo = s[:mid]
    hi = s[mid + (n % 2):]
    q1 = lo[len(lo) // 2]
    q3 = hi[len(hi) // 2]
    return q3 - q1

def _flush_stdin():
    """Discard any keystrokes that accumulated in stdin during long computations."""
    while msvcrt.kbhit():
        msvcrt.getch()

def _norm_hdr(col):
    h = col.lower()
    h = re.sub(r'[/\|\u2022\xb7\u2013\u2014\-]', '_', h)
    h = re.sub(r'[(){}\[\]]', '_', h)
    h = re.sub(r'[\u00b0\'\"\.,;!?]', '', h)
    h = re.sub(r'[\s_]+', '_', h)
    return h.strip('_')
def _has_kw(col, kw_list):
    h = _norm_hdr(col)
    return any(kw in h for kw in kw_list)

def _is_numeric_str(v):
    if v is None: return False
    if isinstance(v, (int, float)): return not (isinstance(v,float) and math.isnan(v))
    s = str(v).strip()
    if not s: return False
    try: float(s); return True
    except: return False

def _is_null_str(v):
    if v is None: return True
    if isinstance(v, float) and math.isnan(v): return True
    NULL = {"","na","n/a","nan","null","none","-","--",".","..","#n/a",
            "#value!","#ref!","#div/0!","#null!","#name?","#num!"}
    return str(v).strip().lower() in NULL

def _looks_like_int_id(vals, max_unique=20):
    """Return True if vals look like a small-cardinality integer identifier
    (tool ID, chamber number, step number, slot, etc.).
    All values must be whole numbers and unique count must be â¤ max_unique.
    """
    nums = [v for v in vals if _is_numeric_str(v)]
    if not nums: return False
    if len(set(str(v) for v in nums)) > max_unique: return False
    for v in nums:
        try:
            fv = float(str(v).strip())
            if abs(fv - round(fv)) > 0.001: return False
        except: return False
    return True


def _looks_like_run_id(vals):
    """Return True if vals look like a run/lot/batch identifier column.
    These are integer IDs that repeat in contiguous blocks (one value per run)
    but can have many more unique values than a tool/step ID.
    Criteria:
      - All non-null values are whole numbers
      - Values appear in contiguous blocks (block_ratio â¤ 2)
      - At least 2 unique values and at least 3 rows per unique value on average
    """
    nums = [v for v in vals if _is_numeric_str(v)]
    if not nums: return False
    for v in nums:
        try:
            fv = float(str(v).strip())
            if abs(fv - round(fv)) > 0.001: return False
        except: return False
    sv = [str(v) for v in nums]
    n_unique = len(set(sv))
    n_vals   = len(sv)
    if n_unique < 2: return False
    if n_vals / n_unique < 3: return False   # need â¥3 rows per run on average
    transitions = sum(1 for i in range(1, n_vals) if sv[i] != sv[i-1])
    block_ratio = transitions / max(n_unique, 1)
    return block_ratio <= 2   # each value appears in â¤2 contiguous blocks

def _score_run_candidate(vals):
    if not vals: return 0
    sv = [str(v) for v in vals]
    n_unique = len(set(sv)); n_vals = len(sv)
    if n_unique < 2: return 0
    if n_unique > n_vals * 0.5: return 0
    transitions = sum(1 for i in range(1,n_vals) if sv[i] != sv[i-1])
    block_ratio = transitions / max(n_unique, 1)
    if block_ratio > 5: return 0
    if n_vals / n_unique < 3: return 0
    return (n_unique ** 2) / (1 + block_ratio)

def discover_schema(rows, headers):
    """
    Auto-discover all structural columns from the data.
    Returns a dict with keys:
      time_col        -- recipe elapsed time (x-axis)
      run_col         -- column that identifies individual runs
      group_cols      -- columns that identify equipment groups (tool, chamberâ¦)
      sensor_cols     -- numeric measurement columns
      step_name_col   -- step name column (or None)
      step_num_col    -- step number column (or None)
      step_elapsed_col-- step-level elapsed time (or None)
      step_seq_col    -- sequence column for step ordering (e.g. CalcStepSeq)
      alarm_code_col  -- alarm code column (or None)
      alarm_status_col-- alarm status/active column (or None)
      event_name_col  -- event name column (or None)
      event_desc_col  -- event description column (or None)
      id_cols         -- all identifier columns
    """
    if not rows:
        return {}

    # ââ per-column stats ââ
    col_info = {}
    for h in headers:
        if not h: continue
        vals     = [r.get(h) for r in rows]
        non_null = [v for v in vals if not _is_null_str(v)]
        nn = len(non_null)
        if nn == 0:
            col_info[h] = dict(numeric_frac=0, n_unique=0, nn=0)
            continue
        num_count = sum(1 for v in non_null if _is_numeric_str(v))
        n_unique  = len(set(str(v) for v in non_null))
        col_info[h] = dict(numeric_frac=num_count/nn, n_unique=n_unique, nn=nn)

    # ââ time column: most monotonically-increasing numeric column
    #    whose name contains a time keyword ââ
    time_col = None
    best_t_score = -1
    for h, info in col_info.items():
        if not _has_kw(h, TIME_KEYWORDS): continue
        if info["numeric_frac"] < 0.5: continue
        float_vals = []
        for r in rows:
            v = r.get(h)
            if _is_numeric_str(v):
                try: float_vals.append(float(str(v).strip()))
                except: pass
        if len(float_vals) < 2: continue
        inc = sum(1 for i in range(1,len(float_vals)) if float_vals[i] > float_vals[i-1])
        flat= sum(1 for i in range(1,len(float_vals)) if abs(float_vals[i]-float_vals[i-1])<1e-9)
        inc_frac  = inc  / max(len(float_vals)-1,1)
        flat_frac = flat / max(len(float_vals)-1,1)
        if flat_frac > 0.7: continue   # constant per run â not elapsed
        score = info["numeric_frac"]*5 + (10 if _has_kw(h,["elapsed"]) else 0) + inc_frac*3
        if score > best_t_score:
            best_t_score = score; time_col = h

    # ââ classify remaining columns ââ
    # This replicates the proven logic from blind_lead_lag_detection_standalone.py
    id_cols, sensor_cols, skip_cols = [], [], []
    time_consumed = {time_col} if time_col else set()

    for h, info in col_info.items():
        if not h or h in time_consumed: continue
        is_num = info["numeric_frac"] >= NUMERIC_FRAC_THRESHOLD
        non_null_vals = [r.get(h) for r in rows if not _is_null_str(r.get(h))]

        # Rule 1: free-text skip columns
        if _has_kw(h, SKIP_KEYWORDS) and not is_num:
            skip_cols.append(h); continue

        # Rule 2: non-numeric â ID
        if not is_num:
            id_cols.append(h); continue

        # ââ numeric column from here ââ

        # Rule 3a: ID keyword + small-cardinality integer (tool, slot, step, etc.)
        if _has_kw(h, ID_KEYWORDS) and _looks_like_int_id(non_null_vals, max_unique=20):
            id_cols.append(h); continue

        # Rule 3b: ID keyword + looks like a run/lot/batch integer ID
        #   (many unique integer values appearing in contiguous blocks,
        #    e.g. Run column with 109 unique run IDs each spanning ~8k rows)
        if _has_kw(h, ID_KEYWORDS) and _looks_like_run_id(non_null_vals):
            id_cols.append(h); continue

        # Rule 4: label keyword + very few unique values (â¤ 10)
        if _has_kw(h, LABEL_KEYWORDS) and info["n_unique"] <= 10:
            id_cols.append(h); continue

        # Rule 5: constant columns (â¤ 2 unique values across entire dataset)
        #   -- these carry no information for anomaly detection
        #   (e.g. UVAErrorCount always 0, HasComments always 0)
        if info["n_unique"] <= 2:
            id_cols.append(h); continue

        # Rule 5b: event-counter / tally columns -- discrete integer tallies
        #   that are not continuous sensor signals.  Excluding these prevents
        #   a cumulative count from dominating the anomaly score and being
        #   mis-labelled as the "lead sensor".
        if _has_kw(h, EVT_COUNT_KW):
            id_cols.append(h); continue

        # Rule 6: DEFAULT -- numeric column with meaningful variation â SENSOR
        sensor_cols.append(h)

    # -- run column --
    # Priority 1: explicit run-ID column names (Run, RunID, Run_ID, Run_Number ...)
    #   Prefer columns whose normalised header exactly matches a known run-ID
    #   keyword AND whose values pass the run/batch shape test.  This prevents
    #   SUBSTRATE_ID (high cardinality, 1 row per wafer) from winning the scorer.
    RUN_EXPLICIT_KW = {'run_id', 'runid', 'run_number', 'runnumber',
                       'run_no', 'runno', 'run_seq', 'run'}
    run_col = None
    best_r  = -1
    for col in (id_cols + sensor_cols):
        if col == time_col: continue
        hn = _norm_hdr(col)
        if hn not in RUN_EXPLICIT_KW: continue
        vals = [r.get(col) for r in rows if r.get(col) is not None]
        if not vals: continue
        if not _looks_like_run_id(vals): continue
        sc = _score_run_candidate(vals)
        if sc > best_r: best_r = sc; run_col = col

    # Priority 2: generic scoring across all id_cols (original logic)
    if run_col is None:
        for col in id_cols:
            if col == time_col: continue
            vals = [r.get(col) for r in rows if r.get(col) is not None]
            sc = _score_run_candidate(vals)
            if sc > best_r: best_r = sc; run_col = col
        # fallback: sensors that look like run/lot/batch integer IDs
        for col in sensor_cols:
            vals = [r.get(col) for r in rows if r.get(col) is not None]
            if not vals: continue
            if not _looks_like_run_id(vals): continue
            sc = _score_run_candidate(vals)
            if sc > best_r: best_r = sc; run_col = col

    # Priority 3: time-based run estimation
    #   No explicit run column found -- detect run boundaries from RunStartTime
    #   transitions (Option A) or from large gaps in the wall-clock timestamp
    #   (Option B, last resort).
    if run_col is None:
        _hdr_lower = {h.lower().replace('_', '').replace(' ', ''): h
                      for h in (list(rows[0].keys()) if rows else [])}
        # Option A: RunStartTime changes signal a new run
        _run_start_cand = None
        for _c in ['runstarttime', 'starttime', 'run_start', 'start_time',
                   'runstart', 'jobstarttime', 'processstart']:
            if _c in _hdr_lower:
                _run_start_cand = _hdr_lower[_c]; break
        if _run_start_cand:
            _lbl = {}; _cnt = 0; _prev = None
            for r in rows:
                rst = str(r.get(_run_start_cand) or '').strip()
                if rst and rst != _prev:
                    if rst not in _lbl:
                        _cnt += 1; _lbl[rst] = _cnt
                    _prev = rst
            if _cnt >= 2:
                _prev = None
                for r in rows:
                    rst = str(r.get(_run_start_cand) or '').strip()
                    if rst: _prev = rst
                    r['_est_run_id'] = str(_lbl.get(_prev or '', 1))
                run_col = '_est_run_id'
                print(f'  No run column found -- estimated {_cnt} runs '
                      f'from {repr(_run_start_cand)} transitions')
        # Option B: large gaps in the wall-clock timestamp column
        if run_col is None:
            _ts_cand = None
            for _c in ['timestamp', 'ts', 'datetime', 'time']:
                if _c in _hdr_lower:
                    _ts_cand = _hdr_lower[_c]; break
            if _ts_cand:
                import datetime as _dtm
                _tv = []
                for r in rows:
                    s = str(r.get(_ts_cand) or '').strip(); _p = None
                    for fmt in ('%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S',
                                '%Y/%m/%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S.%f',
                                '%Y/%m/%d %H:%M:%S', '%Y-%m-%d %H:%M:%S',
                                '%m/%d/%Y %H:%M:%S', '%d/%m/%Y %H:%M:%S'):
                        try: _p = _dtm.datetime.strptime(s, fmt).timestamp(); break
                        except (ValueError, TypeError): continue
                    _tv.append(_p)
                _gaps = sorted(abs(_tv[i] - _tv[i-1])
                               for i in range(1, len(_tv))
                               if _tv[i] is not None and _tv[i-1] is not None
                               and abs(_tv[i] - _tv[i-1]) > 0)
                if _gaps:
                    _thresh = max(60.0, _gaps[len(_gaps)//2] * 20)
                    _rid = 1
                    for i, r in enumerate(rows):
                        if (i > 0 and _tv[i] is not None and _tv[i-1] is not None
                                and (_tv[i] - _tv[i-1]) > _thresh):
                            _rid += 1
                        r['_est_run_id'] = str(_rid)
                    if _rid >= 2:
                        run_col = '_est_run_id'
                        print(f'  No run column found -- estimated {_rid} runs '
                              f'from timestamp gaps (threshold {_thresh:.0f}s)')

    # -- group columns --
    # ToolName-like columns need a DIFFERENT validation path from other id_cols.
    # When the same run_id appears in two tools (e.g. run 107891 in both
    # 5CAP12N28-Carme-DTC_CHA_3 and 5CAPN12N28-Carme-DTC_CHA_4), the standard
    # "constant within run_id group" check FAILS because the merged run group
    # contains rows from both tools.  ToolName is therefore never added to
    # group_cols and tool_col falls back to a wrong column (e.g. UVANormalCount).
    #
    # Tool columns: validate by counting distinct run_ids owned by each tool.
    # Other id_cols: keep the original constant-per-run check.
    TOOL_NAME_KW = ['toolname', 'tool_name', 'equipmentname', 'equipment_name',
                    'eqpname', 'eqp_name', 'chamber_name']
    group_cols = []
    if run_col:
        run_groups = defaultdict(list)
        for r in rows:
            rk = str(r.get(run_col, '?'))
            run_groups[rk].append(r)
        n_runs = len(run_groups)

        # ── Tool-name candidates ─────────────────────────────────────────────
        # Skips constant-per-run check (would fail for cross-tool run_ids).
        # Instead counts distinct run_ids per tool value.
        _tool_cands = [c for c in id_cols
                       if c != run_col
                       and not _has_kw(c, LABEL_KEYWORDS)
                       and _has_kw(c, TOOL_NAME_KW)]
        for col in _tool_cands:
            vals = [r.get(col) for r in rows if r.get(col) is not None]
            n_unique = len(set(str(v) for v in vals))
            if n_unique < 2 or n_unique > 50: continue
            # Count how many distinct run_ids each tool value has
            runs_per_tool = defaultdict(set)
            for r in rows:
                tv = str(r.get(col, '') or '').strip()
                rv = str(r.get(run_col, '') or '').strip()
                if tv and rv:
                    runs_per_tool[tv].add(rv)
            if not runs_per_tool: continue
            min_runs_tool = max(1, n_runs * 0.05)
            if min(len(rids) for rids in runs_per_tool.values()) >= min_runs_tool:
                group_cols.append(col)
                runs_info = ', '.join(f'{t}: {len(r)} runs' for t, r in sorted(runs_per_tool.items()))
                print(f'  Tool column detected: {col!r} ({n_unique} tools -- {runs_info})')

        # ── Other id_col candidates ──────────────────────────────────────────
        # Must be constant within each run_id group (original logic).
        _other_cands = [c for c in id_cols
                        if c != run_col
                        and not _has_kw(c, LABEL_KEYWORDS)
                        and not _has_kw(c, TOOL_NAME_KW)
                        and _has_kw(c, ID_KEYWORDS)]
        for col in _other_cands:
            vals = [r.get(col) for r in rows if r.get(col) is not None]
            n_unique = len(set(str(v) for v in vals))
            if n_unique < 2 or n_unique > 10: continue
            constant = all(
                len(set(str(r.get(col, '') or '') for r in rrows)) <= 1
                for rrows in run_groups.values())
            if not constant: continue
            runs_per_grp = defaultdict(int)
            for rk, rrows in run_groups.items():
                gv = str(rrows[0].get(col, '?'))
                runs_per_grp[gv] += 1
            if min(runs_per_grp.values()) >= max(2, n_runs * 0.15):
                group_cols.append(col)

    # ââ step columns (exact or substring match, first found wins) ââ
    used = set()
    def _find_col(kw_list):
        # exact normalized match first
        for h in headers:
            if h in used: continue
            hn = _norm_hdr(h)
            for kw in kw_list:
                if hn == kw: used.add(h); return h
        # substring match
        for h in headers:
            if h in used: continue
            hn = _norm_hdr(h)
            for kw in kw_list:
                if kw in hn: used.add(h); return h
        return None

    step_num_col    = _find_col(STEP_NUM_KW)
    step_elapsed_col= _find_col(STEP_ELAPSED_KW)
    step_seq_col    = _find_col(STEP_SEQ_KW)
    step_name_col   = _find_col(STEP_NAME_KW)

    # ââ alarm / event columns ââ
    alarm_code_col   = _find_col(["alarmcode", "alarm_code"])
    alarm_status_col = _find_col(["alarmstatus", "alarm_status",
                                   "alarm_active", "alarmactive"])
    event_name_col   = _find_col(["eventname", "event_name", "event"])
    event_desc_col   = _find_col(["eventdescription", "event_description",
                                   "event_desc", "description"])
    event_source_col = _find_col(["eventsource", "event_source", "source"])

    # ââ ARC / process-event-count columns ââ
    # These are excluded from sensor_cols (integer state tallies, not continuous
    # signals), but the first timestep where they transition 0 â >0 is injected
    # as a synthetic alarm event so it appears in the cascade timeline.
    # Restricted to columns whose name contains arc/RF-specific keywords so that
    # general quality metrics (UVAErrorCount, UVAWarningCount, etc.) are not
    # incorrectly treated as process-level arc events.
    ARC_SPECIFIC_KW = [
        "arc", "arcevt", "arc_evt", "ioname_arc",
        "arc_count", "arccount", "arc_event",
        "refl_pwr_fault", "rf_fault",
    ]
    arc_count_cols = [
        h for h in headers
        if _has_kw(h, EVT_COUNT_KW)
        and _has_kw(h, ARC_SPECIFIC_KW)
        and col_info.get(h, {}).get("numeric_frac", 0) >= 0.5
    ]

    # ââ Final sensor_cols cleanup ââ
    # Collect all structural column names discovered above so they are
    # never treated as sensors, regardless of their numeric content.
    structural = ({
        time_col, run_col,
        step_num_col, step_elapsed_col, step_seq_col, step_name_col,
        alarm_code_col, alarm_status_col, event_name_col, event_desc_col,
        event_source_col,
    } | set(group_cols)) - {None}

    # Also explicitly exclude known recipe-bookkeeping keyword patterns.
    # NOTE: "count" alone is intentionally NOT in this list -- arc event
    # counters (e.g. HighFreqRF_rArcEventCount) are genuine process sensors
    # that are non-zero only on bad runs.  Only exclude bookkeeping counters
    # that are identified by their full compound keyword.
    EXCLUDE_KW = [
        "uva", "quality", "ratio", "seq", "loop",
        "has_", "hascomments", "wafer_number", "accum",
        "uvacount", "uvacritical", "uvaerror", "uvanormal", "uvaoutlier",
        "uvawarning", "totalcount", "errorcount", "warningcount",
    ]

    sensor_cols = [
        s for s in sensor_cols
        if s not in structural
        and not s.endswith("_SP")
        and not any(kw in _norm_hdr(s) for kw in EXCLUDE_KW)
    ]

    # ââ Noise filter: remove sensors that zigzag randomly across all runs ââ
    # A noisy sensor alternates up-down-up rapidly (high zigzag ratio) and
    # has meaningful variation (useful_frac > 0) -- i.e. it's not just flat.
    # Such sensors carry no process information and corrupt lead/lag detection.
    # Metric: zigzag_ratio = fraction of consecutive triples (a,b,c) where
    #         b is a local extremum. Threshold: >0.35 across all runs.
    #         useful_frac  = fraction of timesteps NOT at the modal value.
    #         Both conditions required to flag as noisy (avoids flat sensors).
    # Per-step noise filter: for each (sensor, step) pair compute
    #   noise_ratio  = median(within-run std) / step_range
    #   inactive     = abs(median run-mean) / global_max < STEP_INACTIVE_FRAC
    # A (sensor, step) is noisy if noise_ratio > STEP_NOISE_RATIO OR inactive.
    # A sensor is excluded globally if noisy in > STEP_NOISY_STEP_FRAC of steps;
    # otherwise it is only suppressed in the specific steps where it is noisy.
    # The per-step exclusion set is stored in the schema for score_steps.

    _step_col_for_noise = step_seq_col or step_num_col

    run_all_rows_map = defaultdict(list)
    for r in rows:
        rid = str(r.get(run_col, "")) if run_col else "ALL"
        run_all_rows_map[rid].append(r)

    # Global max per sensor (for inactive check)
    sensor_global_max = {}
    for s in sensor_cols:
        all_v = [_f(r.get(s)) for r in rows if _f(r.get(s)) is not None]
        if all_v:
            sensor_global_max[s] = max(abs(v) for v in all_v)

    # Per (sensor, step): collect (within-run std, run_mean, step_range) tuples
    step_sensor_data = defaultdict(lambda: defaultdict(list))
    for rid, rrows_r in run_all_rows_map.items():
        by_step_vals = defaultdict(lambda: defaultdict(list))
        for r in rrows_r:
            sk = str(r.get(_step_col_for_noise, "ALL")).strip() if _step_col_for_noise else "ALL"
            for s in sensor_cols:
                sv = _f(r.get(s))
                if sv is not None:
                    by_step_vals[sk][s].append(sv)
        for sk, s_dict in by_step_vals.items():
            for s, svals in s_dict.items():
                if len(svals) < 5:
                    continue
                mn_s  = sum(svals) / len(svals)
                std_s = math.sqrt(sum((v - mn_s) ** 2 for v in svals) / len(svals))
                rng_s = max(svals) - min(svals)
                # Second-half stats: skip the first half to avoid counting
                # the initial ramp/settling transition as noise.
                half  = max(1, len(svals) // 2)
                sv2   = svals[half:]
                n2    = len(sv2)
                mn_s2 = sum(sv2) / n2
                std_s2= math.sqrt(sum((v - mn_s2) ** 2 for v in sv2) / n2) if n2 > 1 else 0.0
                rng_s2= max(sv2) - min(sv2)
                # Spike fraction: proportion of stable values that are
                # far from the stable median (catches impulse-noise sensors)
                sv2_s = sorted(sv2)
                med_s2 = sv2_s[n2 // 2]
                spike_thresh = max(abs(med_s2) * 3, rng_s2 * 0.30, 1e-6)
                spk_frac = sum(1 for v in sv2 if abs(v - med_s2) > spike_thresh) / n2
                # Tuple: (full_std, full_mean, full_range,
                #         stable_std, stable_mean, stable_range, stable_spike_frac)
                step_sensor_data[sk][s].append(
                    (std_s, mn_s, rng_s, std_s2, mn_s2, rng_s2, spk_frac))

    noisy_sensor_steps = set()
    step_noisy_count   = defaultdict(int)
    step_total_count   = defaultdict(int)
    all_scored_steps   = [sk for sk in step_sensor_data if sk not in ("-1", "0")]

    for sk in all_scored_steps:
        for s, run_tuples in step_sensor_data[sk].items():
            if len(run_tuples) < 3:
                continue
            step_total_count[s] += 1
            g_max = sensor_global_max.get(s, 0.0)

            # Use STABLE (second-half) statistics so sensors that ramp then
            # settle are not falsely classified as noisy.
            stab_stds  = sorted(t[3] for t in run_tuples)
            stab_means = sorted(t[4] for t in run_tuples)
            stab_rngs  = sorted(t[5] for t in run_tuples)
            stab_spks  = sorted(t[6] for t in run_tuples)
            med_stab_std  = stab_stds[len(stab_stds) // 2]
            med_stab_mean = stab_means[len(stab_means) // 2]
            med_stab_rng  = stab_rngs[len(stab_rngs) // 2]
            med_stab_spk  = stab_spks[len(stab_spks) // 2]

            # Criterion A: high RANDOM noise (std/range > threshold) in stable
            # portion, confirmed by non-trivial spike fraction.
            # Pure process ramps (heater going 14%->60%) have high std/range
            # but zero spike fraction -- they are NOT noise and must be kept.
            # Only flag Criterion A when spikes are also present, so genuine
            # ramp steps are never misclassified as noisy.
            _range_floor = max(abs(med_stab_mean) * 1e-3, 1e-6)
            stab_nr = (med_stab_std / med_stab_rng
                       if med_stab_rng > _range_floor else 0.0)
            noisy_A = stab_nr > STEP_NOISE_RATIO and med_stab_spk > 0.01
            # Criterion B: impulse/spike noise alone (even when std/range is low)
            noisy_C = med_stab_spk > 0.05
            # Criterion D: sensor inactive -- stable mean near-zero vs global max
            noisy_B = (g_max > 1e-9 and abs(med_stab_mean) < STEP_INACTIVE_FRAC * g_max)

            if noisy_A or noisy_B or noisy_C:
                noisy_sensor_steps.add((s, sk))
                step_noisy_count[s] += 1

    noisy_sensors = set()
    for s in sensor_cols:
        total = step_total_count.get(s, 0)
        noisy = step_noisy_count.get(s, 0)
        if total >= 3 and noisy / total > STEP_NOISY_STEP_FRAC:
            sp_col_candidate = s + "_SP"
            if sp_col_candidate in (rows[0].keys() if rows else []):
                sp_pcts = []
                for rrows_r in run_all_rows_map.values():
                    for r in rrows_r:
                        va = _f(r.get(s))
                        vs = _f(r.get(sp_col_candidate))
                        if va is not None and vs is not None and vs != 0:
                            sp_pcts.append(abs(va - vs) / abs(vs))
                if sp_pcts:
                    sp_pcts.sort()
                    if sp_pcts[len(sp_pcts) // 2] < 0.05:
                        continue
            # Do not globally exclude a sensor that has at least one step
            # where the STABLE portion passes all three noise checks:
            #   - stable noise ratio (std/range) < STEP_NOISE_RATIO
            #   - spike fraction < 5%
            #   - sensor is not inactive (mean > STEP_INACTIVE_FRAC * global_max)
            has_clean_stable_step = False
            g_max_s = sensor_global_max.get(s, 0.0)
            for sk2 in all_scored_steps:
                if (s, sk2) in noisy_sensor_steps:
                    continue
                tups2 = step_sensor_data[sk2].get(s, [])
                if len(tups2) < 3:
                    continue
                ss2   = sorted(t[3] for t in tups2)
                sm2   = sorted(t[4] for t in tups2)
                sr2   = sorted(t[5] for t in tups2)
                spk2  = sorted(t[6] for t in tups2)
                med_ss2  = ss2[len(ss2) // 2]
                med_sm2  = sm2[len(sm2) // 2]
                med_sr2  = sr2[len(sr2) // 2]
                med_spk2 = spk2[len(spk2) // 2]
                _rf2 = max(abs(med_sm2) * 1e-3, 1e-6)
                nr2  = med_ss2 / med_sr2 if med_sr2 > _rf2 else 0.0
                inact2 = g_max_s > 1e-9 and abs(med_sm2) < STEP_INACTIVE_FRAC * g_max_s
                # Mirror the noisy_A criterion: high std/range is only noise when
                # spikes are also present; ramps pass this check.
                noisy_A2 = nr2 > STEP_NOISE_RATIO and med_spk2 > 0.01
                noisy_C2 = med_spk2 > 0.05
                if (not noisy_A2 and not noisy_C2 and med_sr2 > 1e-6 and not inact2):
                    has_clean_stable_step = True
                    break
            if has_clean_stable_step:
                continue
            noisy_sensors.add(s)

    if noisy_sensors:
        print(f"  Noisy sensors excluded globally ({len(noisy_sensors)}): "
              f"{sorted(noisy_sensors)[:8]}{'...' if len(noisy_sensors) > 8 else ''}")
    n_step_excl = sum(1 for (s, _) in noisy_sensor_steps if s not in noisy_sensors)
    if n_step_excl:
        print(f"  Per-step noise exclusions (sensor noisy in specific steps only): {n_step_excl}")

    sensor_cols = [s for s in sensor_cols if s not in noisy_sensors]

    # ââ Flow sensor filter: inactive gas exclusion ââ
    #
    #    Step 1 -- check only columns whose name contains "flow" (case-insensitive).
    #    Compute the median per-run range (max â min within each run) across all
    #    runs.  If that median is < 10 sccm the gas is considered inactive for
    #    this recipe (always-off or at a fixed setpoint with negligible variation).
    #
    #    Step 2 -- for each inactive gas, extract its "gas token": the substring
    #    of the column name between the MFC/Mfc prefix and the Flow/flow suffix.
    #    Drop ALL sensor columns that share that same MFC-gas prefix -- this
    #    removes the associated Pressure, Temperature, SP, Ratio, etc. columns
    #    for that gas as well.
    #
    #    Example: MFC_AR_Flow  median_range = 0.6 sccm  â  gas token = "MFC_AR"
    #             â also drop MFC_AR_Pressure, MFC_AR_Temperature, MFC_AR_Flow_SP
    #
    #    For non-MFC flow columns (Gas_Total_Flow, O2CLN_Mfc_rFlow, etc.) that
    #    are flat, only that specific column is dropped -- there is no gas-token
    #    cascade.
    #
    #    Columns that contain "flow" but have no MFC/Mfc in their name AND have
    #    median range < 10 sccm are dropped individually (no cascade needed since
    #    there are no sibling MFC columns to clean up).
    FLOW_MIN_RANGE_SCCM = 10.0   # sccm

    def _flow_median_range(col):
        """Return median per-run range for a column, or None if no data."""
        per_run_ranges = []
        for rrows_r in run_all_rows_map.values():
            vals = [_f(r.get(col)) for r in rrows_r if _f(r.get(col)) is not None]
            if len(vals) >= 2:
                per_run_ranges.append(max(vals) - min(vals))
        if not per_run_ranges:
            return None
        per_run_ranges.sort()
        return per_run_ranges[len(per_run_ranges) // 2]

    def _gas_token(col):
        """Extract the gas token from a Flow column name.

        For MFC-style names:  MFC_<GAS_TOKEN>_Flow[_SP]
          e.g. MFC_AR_Flow        â "MFC_AR"
               MFC_AR_CLN_Flow    â "MFC_AR_CLN"
               MFC_NH3_LO_Flow_SP â "MFC_NH3_LO"
               O2CLN_Mfc_rFlow    â "O2CLN_Mfc"

        Only columns whose name starts directly with MFC_ or contains _Mfc_
        are considered primary flow readings.  Columns starting with VS_ or
        other validation/ratio prefixes are excluded -- they are derived metrics,
        not primary MFC flow sensors, and should not drive gas-token cascades.

        Returns None for flow columns with no primary-MFC component.
        """
        # Exclude validation/ratio prefixes -- not primary MFC flow readings
        if re.match(r'VS_', col, re.IGNORECASE):
            return None
        # find the 'flow' word boundary -- also handles rFlow/fFlow MFC naming
        # where a single lowercase r or f precedes "flow" (e.g. O2CLN_Mfc_rFlow)
        m = re.search(r'(?<![a-zA-Z0-9])[rf]?flow(?![a-zA-Z0-9])', col, re.IGNORECASE)
        if m is None:
            return None
        # Strip single leading r/f from the match to get the true prefix boundary
        match_str = m.group(0)
        flow_start = m.start()
        if match_str[0].lower() in ('r', 'f') and len(match_str) > 4:
            flow_start += 1  # skip the r/f prefix character
        prefix = col[:flow_start].rstrip('_')
        if not re.search(r'mfc', prefix, re.IGNORECASE):
            return None
        return prefix   # e.g. "MFC_AR", "MFC_AR_CLN", "O2CLN_Mfc"

    # Pass 1: evaluate every Flow column in sensor_cols.
    # Build a map of gas_token â median_range so we can resolve conflicts
    # where a shorter token (MFC_AR) and a longer token (MFC_AR_CLN) both exist.
    token_range = {}   # gas_token â median flow range
    flat_flow_cols = set()   # non-MFC flow cols dropped individually

    for s in sensor_cols:
        if 'flow' not in s.lower():
            continue
        med = _flow_median_range(s)
        token = _gas_token(s)
        if token:
            # keep the largest range seen for this token (most conservative)
            if med is not None:
                token_range[token] = max(token_range.get(token, 0.0), med)
        else:
            if med is not None and med < FLOW_MIN_RANGE_SCCM:
                flat_flow_cols.add(s)

    # A token is inactive only if its own flow range is < threshold AND no
    # longer more-specific token that starts with this token is active.
    # Example: MFC_AR (range 0.6) is inactive.  MFC_AR_CLN (range 10006) is
    # active.  MFC_AR_CLN_* columns must NOT be cascaded away by MFC_AR.
    # A token is inactive if its own flow range is < threshold.
    # The "active longer token" protection is applied only during the cascade
    # step below, not here -- so MFC_AR is inactive even though MFC_AR_CLN is
    # active.  This ensures MFC_AR_Flow/Pressure/Temperature are dropped while
    # MFC_AR_CLN_Flow/Pressure/Temperature are protected by the cascade guard.
    inactive_gas_tokens = {
        token for token, rng in token_range.items()
        if rng < FLOW_MIN_RANGE_SCCM
    }

    # Pass 2: cascade -- drop all sensor_cols that share an inactive gas prefix.
    # A column belongs to token T if it starts with T + "_" or equals T,
    # AND there is no longer active token that is a more specific match for it.
    active_tokens = {t for t, r in token_range.items() if r >= FLOW_MIN_RANGE_SCCM}

    mfc_gas_excluded = set()
    for s in sensor_cols:
        for token in inactive_gas_tokens:
            if not (s == token or s.lower().startswith(token.lower() + '_')):
                continue
            # Do not exclude if a longer active token is a better match
            has_active_match = any(
                s == at or s.lower().startswith(at.lower() + '_')
                for at in active_tokens
                if len(at) > len(token)
            )
            if not has_active_match:
                mfc_gas_excluded.add(s)
                break

    flow_flat_sensors = flat_flow_cols | mfc_gas_excluded

    if inactive_gas_tokens:
        print(f"  Inactive gas tokens (flow range < {FLOW_MIN_RANGE_SCCM} sccm): "
              f"{sorted(inactive_gas_tokens)}")
    if flow_flat_sensors:
        print(f"  Inactive-gas sensors excluded ({len(flow_flat_sensors)}): "
              f"{sorted(flow_flat_sensors)[:10]}"
              f"{'...' if len(flow_flat_sensors) > 10 else ''}")
    sensor_cols = [s for s in sensor_cols if s not in flow_flat_sensors]

    # ââ Monotonic drift filter: remove sensors with pure gradual drift ââ
    #    (chamber aging, consumable wear) across ALL runs -- not fault step-changes.
    #
    #    Key distinction:
    #      - Pure slow drift:    monotonic trend visible in BOTH the early half
    #                            AND the late half of runs independently.
    #      - Fault step-change:  early half is flat (corr â 0), late half may
    #                            be flat at a different level.
    #      - Combined drift+fault: treated as a sensor of interest -- keep it.
    #
    #    Criteria (ALL three must be true to exclude):
    #      1. |Spearman corr| > 0.85 across the first 40% of runs (strong early trend)
    #      2. |Spearman corr| > 0.85 across the last  40% of runs (trend continues)
    #      3. The range of per-run means over the early half is > 5% of global range
    #         (actual drift, not just numerical noise)
    #    This is much stricter than the old single-window 0.75 threshold, so
    #    fault sensors that happen to drift slightly don't get excluded.
    DRIFT_CORR_THRESH      = 0.85   # both halves must exceed this
    DRIFT_EARLY_FRAC       = 0.40   # first 40% of runs = "early"
    DRIFT_LATE_FRAC        = 0.40   # last  40% of runs = "late"
    DRIFT_RANGE_FRAC       = 0.05   # early-half drift must be â¥5% of global range
    DRIFT_VERYEARLY_FRAC   = 0.20   # first 20% of runs = "very early" window
    DRIFT_VERYEARLY_THRESH = 0.20   # very-early range must be â¥20% of global range
    #   Rationale: genuine warmup/consumable drift produces large variation in
    #   the very first runs (e.g. rReflPwr goes 1.2â0.55 across all good runs).
    #   Fault step-changes are flat in the first N runs and only jump later --
    #   their very-early range is <5% of global range.  Requiring the very-early
    #   window to account for â¥20% of global range prevents fault sensors
    #   (DCBias, ImpedanceR) from being incorrectly excluded even when their
    #   whole-dataset Spearman happens to be high due to the fault step inflating
    #   the late-half correlation.
    MIN_RUNS_FOR_DRIFT = 10

    rough_scoring_excl = set()   # sensors excluded from stable-core rough pass only

    if run_col:
        # Get run order from the data (sorted by first appearance)
        run_order = []
        seen_runs = set()
        for r in rows:
            rid = str(r.get(run_col, "")).strip()
            if rid and rid not in seen_runs:
                run_order.append(rid)
                seen_runs.add(rid)

        n_early      = max(MIN_RUNS_FOR_DRIFT, int(len(run_order) * DRIFT_EARLY_FRAC))
        n_late       = max(MIN_RUNS_FOR_DRIFT, int(len(run_order) * DRIFT_LATE_FRAC))
        n_veryearly  = max(3, int(len(run_order) * DRIFT_VERYEARLY_FRAC))
        early_runs      = run_order[:n_early]
        late_runs       = run_order[-n_late:]
        veryearly_runs  = run_order[:n_veryearly]

        if len(early_runs) >= MIN_RUNS_FOR_DRIFT:
            drifting_sensors = set()
            rough_scoring_excluded_sensors = set()   # local; merged into rough_scoring_excl below

            def _spearman(positions, values):
                n_v = len(positions)
                if n_v < 3:
                    return 0.0
                def _rank(lst):
                    si = sorted(range(n_v), key=lambda i: lst[i])
                    ranks = [0.0] * n_v
                    for rank, idx in enumerate(si):
                        ranks[idx] = float(rank)
                    return ranks
                rp = _rank(positions); rv = _rank(values)
                mrp = sum(rp) / n_v;   mrv = sum(rv) / n_v
                num  = sum((rp[i]-mrp)*(rv[i]-mrv) for i in range(n_v))
                den1 = math.sqrt(sum((rp[i]-mrp)**2 for i in range(n_v)))
                den2 = math.sqrt(sum((rv[i]-mrv)**2 for i in range(n_v)))
                return num / max(den1*den2, 1e-12)

            def _run_means(run_list, s):
                out = []
                for i, rid in enumerate(run_list):
                    rrows_r = run_all_rows_map.get(rid, [])
                    vals = [_f(r.get(s)) for r in rrows_r if _f(r.get(s)) is not None]
                    if vals:
                        out.append((i, sum(vals) / len(vals)))
                return out

            for s in sensor_cols:
                early_means     = _run_means(early_runs,     s)
                late_means      = _run_means(late_runs,      s)
                veryearly_means = _run_means(veryearly_runs, s)

                if len(early_means) < MIN_RUNS_FOR_DRIFT:
                    continue

                # Spearman on early half
                ep = [i for i, _ in early_means]; ev = [m for _, m in early_means]
                corr_early = _spearman(ep, ev)

                # Spearman on late half
                if len(late_means) >= MIN_RUNS_FOR_DRIFT:
                    lp = [i for i, _ in late_means]; lv = [m for _, m in late_means]
                    corr_late = _spearman(lp, lv)
                else:
                    corr_late = corr_early  # not enough late runs -- use early

                # Range checks
                global_vals = []
                for rrows_r in run_all_rows_map.values():
                    global_vals.extend(_f(r.get(s)) for r in rrows_r
                                       if _f(r.get(s)) is not None)
                global_range     = max(global_vals) - min(global_vals) if global_vals else 0
                early_range      = max(ev) - min(ev) if ev else 0
                range_frac       = early_range / max(global_range, 1e-9)

                # Very-early range: if the first 20% of runs already spans â¥20%
                # of the global sensor range, the drift is a warmup/consumable
                # phenomenon (e.g. rReflPwr decaying from 1.2 to 0.55 in the
                # very first runs).  Fault step-changes are flat in early runs
                # and only step later, so their very-early range is <5% of global.
                ve_vals = [m for _, m in veryearly_means]
                veryearly_range  = max(ve_vals) - min(ve_vals) if len(ve_vals) >= 2 else 0
                veryearly_frac   = veryearly_range / max(global_range, 1e-9)

                # Require BOTH halves to show strong monotonic trend AND the
                # early range to be substantial AND the very-early range to be
                # substantial.  The very-early check is the critical guard
                # against accidentally excluding fault sensors (DCBias, ImpedanceR)
                # whose global Spearman is high because the fault step-change
                # inflates the correlation, but whose very-early range is tiny.
                if (abs(corr_early) > DRIFT_CORR_THRESH
                        and abs(corr_late) > DRIFT_CORR_THRESH
                        and range_frac > DRIFT_RANGE_FRAC
                        and veryearly_frac >= DRIFT_VERYEARLY_THRESH):
                    drifting_sensors.add(s)
                    continue   # no need to check rough-pass below

                # Second tier: sensors with strong EARLY-half drift AND large
                # very-early range, but whose LATE-half Spearman is below the
                # full-exclusion threshold (e.g. the fault bloc disrupts the
                # late-half trend).  These sensors should NOT be used for the
                # stable-core rough-scoring pass, because their cross-run drift
                # causes early runs to score artificially high against a
                # mid-run-dominated baseline, incorrectly placing them outside
                # the stable core.  They are KEPT for final anomaly scoring
                # (the clean stable-core baseline handles them correctly there).
                if (abs(corr_early) > DRIFT_CORR_THRESH
                        and veryearly_frac >= DRIFT_VERYEARLY_THRESH
                        and range_frac > DRIFT_RANGE_FRAC):
                    rough_scoring_excluded_sensors.add(s)

            if drifting_sensors:
                print(f"  Slow-drift sensors excluded ({len(drifting_sensors)}): "
                      f"{sorted(drifting_sensors)[:8]}"
                      f"{'...' if len(drifting_sensors) > 8 else ''}")
            rough_scoring_excl.update(rough_scoring_excluded_sensors - drifting_sensors)
            if rough_scoring_excl:
                print(f"  Early-drift sensors excluded from rough scoring only "
                      f"({len(rough_scoring_excl)}): "
                      f"{sorted(rough_scoring_excl)[:8]}")
            sensor_cols = [s for s in sensor_cols if s not in drifting_sensors]

    schema = dict(
        time_col        = time_col,
        run_col         = run_col,
        group_cols      = group_cols,
        sensor_cols     = sensor_cols,
        id_cols         = id_cols,
        step_name_col   = step_name_col,
        step_num_col    = step_num_col,
        step_elapsed_col= step_elapsed_col,
        step_seq_col    = step_seq_col,
        alarm_code_col  = alarm_code_col,
        alarm_status_col= alarm_status_col,
        event_name_col  = event_name_col,
        event_desc_col  = event_desc_col,
        event_source_col= event_source_col,
        arc_count_cols       = arc_count_cols,
        rough_scoring_excl   = rough_scoring_excl,
        noisy_sensor_steps   = noisy_sensor_steps,
    )

    print("\n  -- Schema Discovery --")
    for k, v in schema.items():
        if isinstance(v, list):
            print(f"    {k:20s}: {v[:8]}{'...' if len(v)>8 else ''}")
        else:
            print(f"    {k:20s}: {v}")
    return schema


# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
# 2. DATA LOADING
# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2  --  DATA LOADING & PREPARATION
# ══════════════════════════════════════════════════════════════════════════════

def load_csv(path):
    """Load a CSV file, trying common encodings."""
    for enc in ("utf-8-sig", "utf-8", "latin-1", "cp1252"):
        try:
            with open(path, "r", encoding=enc) as f:
                return list(csv.DictReader(f))
        except (UnicodeDecodeError, UnicodeError):
            continue
    raise ValueError(f"Could not decode file with any common encoding: {path}")


def _f(v):
    """Safe string â float, returns None on failure."""
    if v is None:
        return None
    s = str(v).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
# 2. DATA PREPARATION
# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ

def prepare_runs(rows, schema):
    """
    Group rows into {tool: {run_id: [rows]}} using auto-discovered columns.
    Injects per-row computed fields:
        _elapsed_s      -- recipe elapsed time in seconds
        _step_elapsed_s -- step-level elapsed time in seconds (or None)
        _step_number    -- step number string
    Collects alarm events and step boundaries.
    Returns (tool_runs, step_info, alarm_info).
    """
    time_col         = schema.get("time_col")
    run_col          = schema.get("run_col")
    group_cols       = schema.get("group_cols", [])
    sensor_cols      = schema.get("sensor_cols", [])
    step_num_col     = schema.get("step_num_col")
    step_elapsed_col = schema.get("step_elapsed_col")
    step_seq_col     = schema.get("step_seq_col")
    step_name_col    = schema.get("step_name_col")
    alarm_code_col   = schema.get("alarm_code_col")
    event_name_col   = schema.get("event_name_col")
    event_desc_col   = schema.get("event_desc_col")
    event_source_col = schema.get("event_source_col")

    # Determine tool column: first group_col, else "ALL"
    tool_col = group_cols[0] if group_cols else None

    # ââ Determine time source and scale ââââââââââââââââââââââââââââââââââââââ
    # Priority:
    #   1. time_col (e.g. RECIPE_ELAPSED_TIME) if it has numeric data
    #   2. Wall-clock elapsed from TimeStamp - RunStartTime
    #   3. Row index within each run (last resort)
    time_vals_sample = []
    for r in rows[:500]:
        v = _f(r.get(time_col, "")) if time_col else None
        if v is not None and v > 0:
            time_vals_sample.append(v)

    # Check if time_col is actually populated
    time_col_empty = len(time_vals_sample) == 0

    time_scale = 1.0
    use_wall_clock = False
    timestamp_col  = None
    run_start_col  = None

    if time_col_empty:
        # Try to find a wall-clock timestamp column and a run-start column
        hdr_lower = {h.lower().replace('_','').replace(' ',''): h
                     for h in (rows[0].keys() if rows else [])}
        for cand in ['timestamp', 'ts', 'datetime', 'time']:
            if cand in hdr_lower:
                timestamp_col = hdr_lower[cand]; break
        for cand in ['runstarttime', 'starttime', 'run_start', 'start_time']:
            if cand in hdr_lower:
                run_start_col = hdr_lower[cand]; break

        if timestamp_col or run_start_col:
            use_wall_clock = True
            print(f"  RECIPE_ELAPSED_TIME is empty -- using wall-clock "
                  f"elapsed from '{timestamp_col or run_start_col}'")
        else:
            print("  WARNING: no elapsed time column found -- using row index as time")
    else:
        max_t = max(time_vals_sample)
        if max_t > 1e5:
            time_scale = 1.0 / 1000.0   # ms â s
        elif max_t > 1000:
            gaps = [abs(time_vals_sample[i]-time_vals_sample[i-1])
                    for i in range(1, len(time_vals_sample))
                    if abs(time_vals_sample[i]-time_vals_sample[i-1]) > 0]
            if gaps:
                med_gap = sorted(gaps)[len(gaps)//2]
                if med_gap > 50:
                    time_scale = 1.0 / 1000.0

    # Pre-parse run-start times for wall-clock mode
    def _parse_ts(s):
        """Parse timestamp string to seconds-since-epoch (float). Returns None on failure."""
        if not s or not s.strip():
            return None
        for fmt in ('%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S',
                    '%Y/%m/%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S.%f',
                    '%Y/%m/%d %H:%M:%S',   '%Y-%m-%d %H:%M:%S',
                    '%m/%d/%Y %H:%M:%S',   '%d/%m/%Y %H:%M:%S'):
            try:
                import time as _time
                import datetime as _dt
                return _dt.datetime.strptime(s.strip(), fmt).timestamp()
            except ValueError:
                continue
        return None

    run_start_ts = {}   # run_id -> start timestamp in seconds
    if use_wall_clock:
        col_for_start = run_start_col or timestamp_col
        for r in rows:
            rid = str(r.get(run_col, "RUN_ALL")).strip() if run_col else "RUN_ALL"
            if rid not in run_start_ts:
                ts = _parse_ts(r.get(col_for_start, ""))
                if ts is not None:
                    run_start_ts[rid] = ts

    # ââ Pre-compute step-start timestamps per (run, step_seq) when no
    #    step_elapsed_col exists.  This lets us derive _step_elapsed_s from
    #    wall-clock timestamps even when the dataset has no explicit column for it.
    #    Strategy: first pass over all rows grouped by run; within each run track
    #    when CalcStepSeq (or step_num_col) changes -- the first wall-clock
    #    timestamp for that step value is the step start time.
    # ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
    step_start_ts_map = {}   # (run_id, step_seq_val) -> wall-clock seconds

    need_step_elapsed_from_ts = (step_elapsed_col is None) and (use_wall_clock or not time_col_empty)
    _step_key_col = step_seq_col or step_num_col   # column that identifies the step

    if need_step_elapsed_from_ts and _step_key_col:
        # Collect first wall-clock ts per (run, step_key) in one pass
        col_for_ts = timestamp_col if use_wall_clock else None
        if col_for_ts is None:
            # Even if we use RECIPE_ELAPSED_TIME for run elapsed, we still need
            # wall-clock timestamps to compute step-relative elapsed.
            hdr_lower2 = {h.lower().replace('_','').replace(' ',''): h
                          for h in (rows[0].keys() if rows else [])}
            for cand in ['timestamp', 'ts', 'datetime', 'time']:
                if cand in hdr_lower2:
                    col_for_ts = hdr_lower2[cand]; break

        if col_for_ts:
            for r in rows:
                rid = str(r.get(run_col, "RUN_ALL")).strip() if run_col else "RUN_ALL"
                sk  = str(r.get(_step_key_col, "")).strip()
                key = (rid, sk)
                if key not in step_start_ts_map:
                    ts_v = _parse_ts(r.get(col_for_ts, ""))
                    if ts_v is not None:
                        step_start_ts_map[key] = ts_v
        else:
            need_step_elapsed_from_ts = False

    # Pre-resolve the wall-clock column for step-elapsed computation once,
    # so we don't repeat the header lookup inside the hot per-row loop.
    _ts_col_for_step_elapsed = None
    if need_step_elapsed_from_ts:
        _ts_col_for_step_elapsed = timestamp_col if use_wall_clock else None
        if _ts_col_for_step_elapsed is None:
            hdr_lower_se = {h.lower().replace('_','').replace(' ',''): h
                            for h in (rows[0].keys() if rows else [])}
            for _cand in ['timestamp', 'ts', 'datetime', 'time']:
                if _cand in hdr_lower_se:
                    _ts_col_for_step_elapsed = hdr_lower_se[_cand]; break
        if _ts_col_for_step_elapsed is None:
            need_step_elapsed_from_ts = False   # no usable timestamp column

    raw       = defaultdict(lambda: defaultdict(list))
    alarm_raw = defaultdict(list)
    run_row_idx = defaultdict(int)   # for row-index fallback

    for r in rows:
        # run identity
        run = str(r.get(run_col, "")).strip() if run_col else ""
        if not run:
            run = "RUN_ALL"
        # tool/group identity
        tool = str(r.get(tool_col, "")).strip() if tool_col else "ALL"
        if not tool:
            tool = "ALL"

        # elapsed time -- three strategies
        if not time_col_empty:
            t_raw = _f(r.get(time_col, "")) if time_col else None
            if t_raw is None or t_raw < 0:
                continue
            elapsed_s = t_raw * time_scale
        elif use_wall_clock:
            ts_str = r.get(timestamp_col, "") if timestamp_col else ""
            ts_abs = _parse_ts(ts_str)
            if ts_abs is None:
                continue
            t_start = run_start_ts.get(run, ts_abs)
            elapsed_s = max(0.0, ts_abs - t_start)
        else:
            # row index within run
            elapsed_s = float(run_row_idx[run])

        run_row_idx[run] += 1
        r["_elapsed_s"] = elapsed_s

        # step fields
        rset = _f(r.get(step_elapsed_col, "")) if step_elapsed_col else None
        if rset is not None:
            r["_step_elapsed_s"] = rset * time_scale
        elif need_step_elapsed_from_ts and _step_key_col:
            # Derive step elapsed from wall-clock timestamps
            sk = str(r.get(_step_key_col, "")).strip()
            step_t0 = step_start_ts_map.get((run, sk))
            if step_t0 is not None and _ts_col_for_step_elapsed:
                ts_v = _parse_ts(r.get(_ts_col_for_step_elapsed, ""))
                r["_step_elapsed_s"] = max(0.0, ts_v - step_t0) if ts_v is not None else None
            else:
                r["_step_elapsed_s"] = None
        else:
            r["_step_elapsed_s"] = None
        r["_step_number"] = str(r.get(step_num_col, "")).strip() if step_num_col else ""

        # alarm/event collection
        alarm_code   = str(r.get(alarm_code_col,   "")).strip() if alarm_code_col   else ""
        event_name   = str(r.get(event_name_col,   "")).strip() if event_name_col   else ""
        event_desc   = str(r.get(event_desc_col,   "")).strip() if event_desc_col   else ""
        event_source = str(r.get(event_source_col, "")).strip() if event_source_col else ""
        step_name    = str(r.get(step_name_col,    "")).strip() if step_name_col    else ""
        is_alarm   = (alarm_code.lower() in ("fault","error","critical","alarm") or
                      "TIMEOUT" in event_name.upper() or
                      "FAILED"  in event_name.upper() or
                      "FAULT"   in event_name.upper())
        if is_alarm and elapsed_s > 0:
            alarm_raw[run].append({
                "elapsed_s":    elapsed_s,
                "alarm_code":   alarm_code,
                "event_name":   event_name,
                "event_source": event_source,
                "description":  event_desc,
                "step_name":    step_name,
            })

        # filter: skip junk rows before the recipe starts
        # When using RECIPE_ELAPSED_TIME: skip t=0 rows with no real step number
        # When using wall-clock/row-index: only skip rows with step number < 1
        sn_int = _safe_int(r["_step_number"], -99)
        if not time_col_empty and elapsed_s <= 0 and sn_int < 1:
            continue
        elif time_col_empty and step_num_col and sn_int < 1:
            continue

        # keep only rows that have at least one numeric sensor value
        if not any(_f(r.get(s)) is not None for s in sensor_cols):
            continue

        raw[tool][run].append(r)

    # sort each run by elapsed time
    tool_runs = {}
    for tool, rdict in raw.items():
        tool_runs[tool] = {}
        for rid, rrows in rdict.items():
            rrows.sort(key=lambda r: r["_elapsed_s"])
            tool_runs[tool][rid] = rrows

    # deduplicate alarms
    alarm_info = {}
    for run_id, alist in alarm_raw.items():
        seen, deduped = set(), []
        for a in sorted(alist, key=lambda x: x["elapsed_s"]):
            key = (round(a["elapsed_s"], 1), a["event_name"])
            if key not in seen:
                seen.add(key); deduped.append(a)
        alarm_info[run_id] = deduped

    # ââ ARC count onset injection ââ
    # For each ARC/event-count column, scan every run's rows in elapsed-time
    # order and inject a synthetic alarm at:
    #   1. The first timestep where the value transitions 0 â >0  ("ARC onset")
    #   2. Each subsequent escalation to a higher integer level (e.g. 1 â 2)
    # These synthetic events appear in the lead/lag cascade alarm timeline and
    # on the SVG charts as vertical markers, telling the analyst exactly when
    # arcing activity began relative to the sensor divergence time.
    arc_count_cols = schema.get("arc_count_cols", [])
    if arc_count_cols:
        for tool, rdict in tool_runs.items():
            for rid, rrows in rdict.items():
                for arc_col in arc_count_cols:
                    prev_level = 0
                    for r in rrows:          # already sorted by elapsed_s
                        v = _f(r.get(arc_col))
                        if v is None:
                            continue
                        level = int(round(v))
                        if level > prev_level:
                            tag = ("ARC_ONSET" if prev_level == 0
                                   else f"ARC_ESCALATE_{prev_level}_TO_{level}")
                            desc = (f"{arc_col} transitioned {prev_level} â {level} "
                                    f"at step {r.get('_step_number','?')}")
                            event = {
                                "elapsed_s":    r["_elapsed_s"],
                                "alarm_code":   "ARC_EVENT",
                                "event_name":   f"ARC.{arc_col}.{tag}",
                                "event_source": arc_col,
                                "description":  desc,
                                "step_name":    (str(r.get(step_name_col, "")).strip()
                                                 if step_name_col else ""),
                                "arc_level":    level,
                                "arc_col":      arc_col,
                                "is_arc_onset": (prev_level == 0),
                            }
                            alarm_info.setdefault(rid, []).append(event)
                            prev_level = level

    # step boundaries using step_seq_col (most reliable) â step_num_col â step_name_col
    step_info = {}
    boundary_col = step_seq_col or step_num_col or step_name_col
    for tool, rdict in tool_runs.items():
        for rid, rrows in rdict.items():
            bounds, prev_val = [], None
            for i, r in enumerate(rrows):
                val = str(r.get(boundary_col, "")).strip() if boundary_col else ""
                if val == prev_val:
                    continue
                bounds.append({
                    "elapsed_s":   r["_elapsed_s"],
                    "step_number": r.get(step_num_col,  "").strip() if step_num_col  else "",
                    "step_name":   r.get(step_name_col, "").strip() if step_name_col else "",
                    "row_idx":     i,
                })
                prev_val = val
            step_info[rid] = bounds

    return tool_runs, step_info, alarm_info


def detect_looped_steps(step_info):
    """{ run_id: set_of_step_names_appearing_more_than_once }"""
    looped = {}
    for rid, bounds in step_info.items():
        cnt = defaultdict(int)
        for b in bounds:
            if b["step_name"]:
                cnt[b["step_name"]] += 1
        looped[rid] = {n for n, c in cnt.items() if c > 1}
    return looped


# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
# 2b. SETPOINT PAIR DETECTION & TWO-FACTOR SENSOR CLASSIFICATION
#
# Two questions per sensor for an anomalous run:
#
#   Factor 1 -- FAR FROM BASELINE?
#     Is the sensor's mean value across this run significantly different from
#     the cross-run baseline (median of all other runs)?
#     Threshold: |run_mean - baseline_median| / baseline_mad > BASELINE_Z_THRESH
#
#   Factor 2 -- FAR FROM ITS OWN SETPOINT (SP)?
#     Is the actual value close to its own SP for this run?
#     Threshold: mean |actual - SP| / |SP| < SP_TRACKING_THRESHOLD
#                i.e. the sensor is tracking its SP well (within 5%)
#
#   Classification:
#     Factor 1 = YES, Factor 2 = YES (close to SP)  â "Change in SP"
#       The run's SP was set to a different value than the other runs.
#       The sensor is following its SP perfectly -- not a fault.
#       Example: Pedestal_Heater_Temperature set to 150Â°C vs 180Â°C in other runs.
#
#     Factor 1 = YES, Factor 2 = NO  (far from SP)  â "Genuine Fault"
#       The sensor deviated from its SP AND from other runs -- real problem.
#       Example: Chamber_Pressure failed to reach LP Clean setpoint.
#
#     Factor 1 = YES, no SP column                  â "Anomalous (no SP)"
#       Sensor differs from baseline but no SP to compare against.
#       Example: Throttle_Valve_Position.
#
#     Factor 1 = NO                                 â "Normal" -- not flagged.
# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ

BASELINE_Z_THRESH    = 1.5    # run mean must be this many robust-sigma from baseline
SP_TRACKING_THRESHOLD = 0.05  # sensor is "tracking SP" if mean pct dev < 5%
SP_OFF_FRAC_THRESHOLD = 0.05  # used for chart shading: off-SP if >5% deviation


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3  --  SP PAIR DETECTION & SENSOR CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════════════

def detect_sp_pairs(all_rows):
    """
    Auto-detect (actual, setpoint) column pairs from CSV headers.

    Supported naming conventions (case-insensitive):
      Suffix  : X_SP, X_Setpoint, X_SetPoint, X_SP_Value, X_Target, X_Desired
      Prefix  : SP_X, Setpoint_X, SetPoint_X, Target_X
      Infix   : X_rSetpoint_Y  (e.g. Chamber_rSetpoint_Temp vs Chamber_rTemp)
      'w' vs 'r' convention: sensor named X_rTemp has SP named X_wTemp
                             (write = commanded, read = actual)

    A pair is confirmed only when both columns have at least some non-zero
    numeric values in the first 200 rows.
    Returns list of (actual_col, sp_col) tuples, sorted by actual_col name.
    """
    if not all_rows:
        return []
    headers = list(all_rows[0].keys())
    header_set = set(headers)

    # Normalise a header to lowercase-underscore for fuzzy matching
    def _n(h):
        import re
        return re.sub(r'[\s\-]+', '_', h.lower())

    # Build a map from normalised name â original header
    norm_to_orig = {_n(h): h for h in headers}

    def _has_data(col):
        return any(_f(r.get(col)) is not None for r in all_rows[:200])

    pairs = []
    seen  = set()   # avoid duplicate (actual, sp) pairs

    for h in headers:
        hn = _n(h)

        # ââ Strategy 1: known SP suffixes appended to actual name ââ
        SP_SUFFIXES = [
            "_sp", "_setpoint", "_set_point", "_sp_value",
            "_target", "_desired", "_cmd", "_command",
            "_setp", "_set_pt", "_spval",
        ]
        for sfx in SP_SUFFIXES:
            if hn.endswith(sfx):
                actual_norm = hn[: -len(sfx)]
                if actual_norm in norm_to_orig:
                    actual_col = norm_to_orig[actual_norm]
                    sp_col     = h
                    key = (actual_col, sp_col)
                    if key not in seen and _has_data(actual_col) and _has_data(sp_col):
                        pairs.append(key)
                        seen.add(key)

        # ââ Strategy 2: known SP prefixes prepended to actual name ââ
        SP_PREFIXES = [
            "sp_", "setpoint_", "set_point_", "target_", "desired_",
            "cmd_", "command_", "setp_",
        ]
        for pfx in SP_PREFIXES:
            if hn.startswith(pfx):
                actual_norm = hn[len(pfx):]
                if actual_norm in norm_to_orig:
                    actual_col = norm_to_orig[actual_norm]
                    sp_col     = h
                    key = (actual_col, sp_col)
                    if key not in seen and _has_data(actual_col) and _has_data(sp_col):
                        pairs.append(key)
                        seen.add(key)

        # ââ Strategy 3: 'w' (write/command) vs 'r' (read/actual) convention ââ
        # e.g. Chamber_wTemp  â  Chamber_rTemp
        #      HighFreqRF_wPower â HighFreqRF_rForwPwr
        # Match: if header contains '_w' followed by a word, look for '_r' + same word
        import re as _re
        m = _re.search(r'(_[wW])([A-Z][A-Za-z0-9]*)', h)
        if m:
            actual_candidate = h[:m.start()] + '_r' + m.group(2)
            if actual_candidate in header_set:
                actual_col = actual_candidate
                sp_col     = h
                key = (actual_col, sp_col)
                if key not in seen and _has_data(actual_col) and _has_data(sp_col):
                    pairs.append(key)
                    seen.add(key)

        # ââ Strategy 4: infix 'Setpoint' / 'SetPoint' in column name ââ
        # e.g. Chamber_SetpointTemp  paired with  Chamber_rTemp
        # Look for any header that contains 'setpoint' and try stripping it
        for sp_kw in ['setpoint', 'set_point', 'setp']:
            if sp_kw in hn:
                actual_norm = hn.replace(sp_kw, '').strip('_')
                if actual_norm in norm_to_orig:
                    actual_col = norm_to_orig[actual_norm]
                    sp_col     = h
                    key = (actual_col, sp_col)
                    if key not in seen and _has_data(actual_col) and _has_data(sp_col):
                        pairs.append(key)
                        seen.add(key)

    return sorted(pairs)


def classify_sensor_deviations(run_dict, sp_pairs, anom_set, sensor_list=None):
    """
    For each anomalous run, classify every sensor (from sensor_list plus
    all SP-pair actuals) into one of:
        "Genuine Fault"   -- far from baseline AND far from its SP
        "Change in SP"    -- far from baseline BUT tracking its SP well
        "Anomalous (no SP)" -- far from baseline, no SP column to compare
        "Normal"          -- not significantly different from baseline

    Returns:
        sensor_classes : { run_id: { sensor: { "category", "bl_z",
                                               "sp_pct_dev", "sp_col" } } }
        sp_stats       : { run_id: { (actual, sp_col): {stats dict} } }
    """
    import math

    # Build a simple cross-run baseline: mean of per-run means (excluding the
    # anomalous run) for each sensor, measured over all timesteps.
    # Use mean-of-means (one value per run) so that run length doesn't bias.
    all_sensors = list(sensor_list) if sensor_list else []
    # add any SP-pair actuals not already in the sensor list
    sp_actuals  = {actual for actual, _ in sp_pairs}
    sp_map      = {actual: sp_col for actual, sp_col in sp_pairs}
    for s in sp_actuals:
        if s not in all_sensors:
            all_sensors.append(s)

    # Per-run mean for each sensor
    run_means = {}   # { rid: { sensor: mean_value } }
    for rid, rrows in run_dict.items():
        run_means[rid] = {}
        for sensor in all_sensors:
            vals = [_f(r.get(sensor)) for r in rrows
                    if _f(r.get(sensor)) is not None
                    and _f(r.get("_elapsed_s", 0)) > 1.0]
            if vals:
                run_means[rid][sensor] = _mean(vals)

    # Compute baseline median + MAD across normal runs
    normal_ids = [rid for rid in run_dict if rid not in anom_set]
    bl_stats   = {}   # { sensor: (median, mad) }
    for sensor in all_sensors:
        norm_means = [run_means[rid][sensor]
                      for rid in normal_ids
                      if sensor in run_means.get(rid, {})]
        if not norm_means:
            continue
        med = _median(norm_means)
        mad = _median([abs(v - med) for v in norm_means]) * 1.4826
        mad = max(mad, abs(med) * 0.01, 1e-6)   # floor
        bl_stats[sensor] = (med, mad)

    # Compute SP tracking stats for every run
    sp_stats = {}
    for rid, rrows in run_dict.items():
        sp_stats[rid] = {}
        for actual, sp_col in sp_pairs:
            active, pct_devs = [], []
            off_ts, off_vs, sp_ts, sp_vs, act_ts, act_vs = [], [], [], [], [], []
            for r in rrows:
                v_act = _f(r.get(actual))
                v_sp  = _f(r.get(sp_col))
                if v_act is None or v_sp is None or v_sp == 0:
                    continue
                t   = r["_elapsed_s"]
                dev = abs(v_act - v_sp)
                pct = dev / abs(v_sp)
                active.append((t, v_act, v_sp))
                pct_devs.append(pct)
                sp_ts.append(t);   sp_vs.append(v_sp)
                act_ts.append(t);  act_vs.append(v_act)
                if pct > SP_OFF_FRAC_THRESHOLD:
                    off_ts.append(t);  off_vs.append(v_act)
            if not active:
                continue
            sp_stats[rid][(actual, sp_col)] = {
                "n_active":  len(active),
                "mean_pct":  _mean(pct_devs),
                "max_pct":   max(pct_devs),
                "frac_off":  len(off_ts) / len(active),
                "mean_dev":  _mean([abs(a - s) for _, a, s in active]),
                "max_dev":   max(abs(a - s) for _, a, s in active),
                "off_ts": off_ts, "off_vs": off_vs,
                "sp_ts":  sp_ts,  "sp_vs":  sp_vs,
                "act_ts": act_ts, "act_vs": act_vs,
            }

    # Classify each sensor for each anomalous run
    sensor_classes = {}
    for rid in anom_set:
        sensor_classes[rid] = {}
        for sensor in all_sensors:
            if sensor not in run_means.get(rid, {}):
                continue
            if sensor not in bl_stats:
                continue
            bl_med, bl_mad = bl_stats[sensor]
            run_mean       = run_means[rid][sensor]
            bl_z           = abs(run_mean - bl_med) / bl_mad

            sp_col     = sp_map.get(sensor)
            sp_pct_dev = None
            if sp_col:
                sd = sp_stats.get(rid, {}).get((sensor, sp_col), {})
                sp_pct_dev = sd.get("mean_pct")   # None if no SP data

            # Classify
            # SP-tracking check takes priority over the z-score check:
            # if the sensor has a setpoint and is following it closely
            # (<5% deviation) it is ALWAYS "Change in SP" -- even if the
            # z-score is low (which happens when the baseline itself is a
            # mix of runs at different setpoints, making the median sit
            # between them and the z-score look artificially small).
            if sp_pct_dev is not None and sp_pct_dev < SP_TRACKING_THRESHOLD:
                category = "Change in SP"
            elif bl_z < BASELINE_Z_THRESH:
                category = "Normal"
            elif sp_pct_dev is None:
                category = "Anomalous (no SP)"
            else:
                category = "Genuine Fault"

            sensor_classes[rid][sensor] = {
                "category":  category,
                "bl_z":      bl_z,
                "run_mean":  run_mean,
                "bl_median": bl_med,
                "sp_col":    sp_col,
                "sp_pct_dev": sp_pct_dev,
            }

    return sensor_classes, sp_stats


def sp_deviation_score(stats):
    """Single score: how far off SP -- used for sorting flagged pairs."""
    return stats.get("frac_off", 0) * stats.get("mean_pct", 0)


# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
# 3. BASELINE  (proven logic from blind_lead_lag_detection_standalone.py)
# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ

def _median(vals):
    if not vals:
        return 0.0
    s = sorted(vals)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n//2-1] + s[n//2]) / 2.0

def _mean(vals):
    return sum(vals) / len(vals) if vals else 0.0

def _safe_int(v, default=-1):
    try: return int(v)
    except: return default

def _nearest(sorted_list, v):
    """Binary-search nearest index.  Works for both float lists and tuple lists."""
    idx = bisect.bisect_left(sorted_list, v)
    return min(idx, len(sorted_list) - 1)


def _row_key(r, BIN=0.5):
    """Return the (step_number, step_elapsed_bin) lookup key for a row."""
    se = r.get("_step_elapsed_s")
    sn = r.get("_step_number", "")
    if se is None or sn == "":
        return None
    return (sn, round(se / BIN) * BIN)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4  --  BASELINE & ANOMALY DETECTION
# ══════════════════════════════════════════════════════════════════════════════

def build_baseline(run_dict, sensor, exclude=None):
    """
    Per-step-timestep median + MAD baseline keyed by (step_number, step_elapsed_bin).

    Returns a dict:
      {
        step_num_str: {
          "bins":    [float, ...],   # sorted step_elapsed bins
          "medians": [float, ...],
          "mads":    [float, ...],
        },
        ...
      }

    Lookup: for a row with (_step_number=sn, _step_elapsed_s=se), find the
    entry for sn and bisect into "bins" -- this guarantees we never cross
    step boundaries when looking up the nearest baseline value.

    Also returns the legacy flat (keys, medians, mads) tuple as a second
    element for backward-compat with z_trace / find_divergence_time callers.
    """
    BIN = 0.5   # seconds
    # by_step[sn][bin] = [values]
    by_step = defaultdict(lambda: defaultdict(list))
    for rid, rrows in run_dict.items():
        if exclude and rid in exclude:
            continue
        for r in rrows:
            v  = _f(r.get(sensor))
            se = r.get("_step_elapsed_s")
            sn = r.get("_step_number", "")
            if v is None or se is None or sn == "":
                continue
            b = round(se / BIN) * BIN
            by_step[sn][b].append(v)

    if not by_step:
        return {}, [], [], []

    step_index = {}   # sn -> {bins, medians, mads}
    all_pos_mads = []

    for sn, bin_vals in by_step.items():
        bins_sorted = sorted(bin_vals.keys())
        step_medians, step_mads = [], []
        for b in bins_sorted:
            vals = bin_vals[b]
            med  = _median(vals)
            mad  = (_median([abs(v - med) for v in vals]) * 1.4826
                    if len(vals) > 1 else 0.0)
            step_medians.append(med)
            step_mads.append(mad)
            if mad > 0:
                all_pos_mads.append(mad)
        step_index[sn] = dict(bins=bins_sorted,
                               medians=step_medians,
                               mads=step_mads)

    global_mad = _median(all_pos_mads) if all_pos_mads else 1.0

    # Apply MAD floor per entry (prevents near-zero MAD from inflating z-scores)
    for sn, entry in step_index.items():
        for i in range(len(entry["mads"])):
            floor = max(abs(entry["medians"][i]) * 0.01,
                        global_mad * 0.1, 1e-9)
            entry["mads"][i] = max(entry["mads"][i], floor)

    # Build legacy flat (keys, medians, mads) for callers that still use it
    def _key_sort(k):
        sn, st = k
        try:    return (int(sn), st)
        except: return (9999,    st)

    flat_keys, flat_meds, flat_mads = [], [], []
    for sn, entry in step_index.items():
        for b, med, mad in zip(entry["bins"], entry["medians"], entry["mads"]):
            flat_keys.append((sn, b))
            flat_meds.append(med)
            flat_mads.append(mad)

    order = sorted(range(len(flat_keys)), key=lambda i: _key_sort(flat_keys[i]))
    flat_keys  = [flat_keys[i]  for i in order]
    flat_meds  = [flat_meds[i]  for i in order]
    flat_mads  = [flat_mads[i]  for i in order]

    return step_index, flat_keys, flat_meds, flat_mads


# -- MFC ramp interval helpers -----------------------------------------------

def _is_mfc_flow_sensor(sensor_name):
    """Return True if this is an MFC flow sensor (ramps should be ignored)."""
    sn = sensor_name.lower()
    return 'mfc' in sn and 'flow' in sn


def _find_sp_col(sensor, row_keys):
    """Return the SP column name for an MFC flow sensor, or None.

    Tries common conventions:
      MFC_AR_CARRIER_Flow  -> MFC_AR_CARRIER_Flow_SP
      MFC_AR_CARRIER_rFlow -> MFC_AR_CARRIER_Flow_SP  (strip leading r/f)
      MFC_AR_CARRIER_fFlow -> MFC_AR_CARRIER_Flow_SP
    Falls back to any column in row_keys whose name contains the same
    gas-token prefix and both 'flow' and 'sp'.
    """
    key_set = set(row_keys)
    # Direct _SP suffix
    if sensor + "_SP" in key_set:
        return sensor + "_SP"
    # Strip leading r/f from the flow suffix (rFlow -> Flow_SP)
    import re
    m = re.search(r'[rf]flow', sensor, re.IGNORECASE)
    if m:
        base = sensor[:m.start()] + "Flow"
        if base + "_SP" in key_set:
            return base + "_SP"
    # Fuzzy: same prefix, contains 'flow' and 'sp'
    sn_lower = sensor.lower()
    prefix = sn_lower[:max(sn_lower.rfind('flow'), 0)].rstrip('_r f'.split()[0])
    for k in key_set:
        kl = k.lower()
        if prefix and kl.startswith(prefix) and 'flow' in kl and 'sp' in kl:
            return k
    return None


def _build_mfc_ramp_intervals(run_dict, sensor, bin_s=1.0):
    """
    Return a list of (t_start, t_end, sp_col, sp_target) tuples where the
    cross-run median of this MFC flow sensor is changing rapidly (a normal ramp).

    sp_col     : name of the setpoint column for this flow sensor (or None)
    sp_target  : median SP value just AFTER the ramp ends (the destination target)

    During scoring, a row inside a ramp window is suppressed only when
    the actual flow is within MFC_OVERSHOOT_FRAC of sp_target.  If the
    flow significantly exceeds sp_target it is still scored as anomalous.
    """
    if not _is_mfc_flow_sensor(sensor):
        return []

    # Find the SP column from any row in run_dict
    sp_col = None
    for rrows in run_dict.values():
        if rrows:
            sp_col = _find_sp_col(sensor, rrows[0].keys())
            break

    # Bin cross-run medians into 1-second buckets
    bin_vals    = defaultdict(list)
    bin_sp_vals = defaultdict(list)
    for rrows in run_dict.values():
        for r in rrows:
            t = r.get("_elapsed_s")
            v = _f(r.get(sensor))
            if t is not None and v is not None:
                b = round(t / bin_s) * bin_s
                bin_vals[b].append(v)
            if sp_col:
                vs = _f(r.get(sp_col))
                if t is not None and vs is not None:
                    bin_sp_vals[round(t / bin_s) * bin_s].append(vs)

    if not bin_vals:
        return []

    times   = sorted(bin_vals)
    medians = [_median(bin_vals[t]) for t in times]

    # Detect bins where the cross-run median changes faster than the threshold
    ramp_bins = set()
    for i in range(1, len(times)):
        dt = times[i] - times[i - 1]
        if dt <= 0:
            continue
        rate = abs(medians[i] - medians[i - 1]) / dt
        if rate > MFC_RAMP_RATE_SCCM_S:
            ramp_bins.add(times[i - 1])
            ramp_bins.add(times[i])

    if not ramp_bins:
        return []

    # Expand ramp bins with padding and merge overlaps
    raw = [(bt - MFC_RAMP_WINDOW_S, bt + MFC_RAMP_WINDOW_S)
           for bt in sorted(ramp_bins)]
    merged = []
    for start, end in sorted(raw):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])

    # For each merged interval compute the SP target (the destination the ramp
    # is heading towards).  Strategy:
    #   1. SP column post-ramp: median SP in a 5-second window after interval
    #   2. SP column during ramp: median SP inside the interval (SP changes
    #      immediately to the new target even while flow is still ramping)
    #   3. Fallback: median flow after the ramp ends
    # The sp_target drives the overshoot check: flow > sp_target*(1+tol) = anomaly.
    result = []
    POST_WIN = 5.0
    for start, end in merged:
        post_bins   = [t for t in times if end          < t <= end + POST_WIN]
        during_bins = [t for t in times if start        <= t <= end           ]

        sp_target = None
        if sp_col:
            # Try post-ramp SP first
            sp_post = [v for b in post_bins   for v in bin_sp_vals.get(b, [])]
            sp_target = _median(sp_post) if sp_post else None
            # Try during-ramp SP if post-ramp is missing or zero
            if sp_target is None or sp_target == 0.0:
                sp_dur = [v for b in during_bins for v in bin_sp_vals.get(b, [])]
                sp_dur_med = _median(sp_dur) if sp_dur else None
                if sp_dur_med is not None and sp_dur_med != 0.0:
                    sp_target = sp_dur_med

        if sp_target is None or sp_target == 0.0:
            # Flow-based fallback: median flow after the ramp
            flow_post = [v for b in post_bins for v in bin_vals.get(b, [])]
            fp_med    = _median(flow_post) if flow_post else None
            if fp_med is not None:
                sp_target = fp_med   # may be 0 for ramp-to-zero -- that's correct

        result.append((start, end, sp_col, sp_target))

    return result


def _in_ramp_interval(t, intervals):
    """Return True if t falls within any ramp interval AND the caller should
    suppress this point (i.e. it is a plain ramp with no overshoot).

    intervals items are (t_start, t_end, sp_col, sp_target) -- the sp fields
    are used by _ramp_should_suppress; this function just checks membership.
    """
    for entry in intervals:
        t_start, t_end = entry[0], entry[1]
        if t_start <= t <= t_end:
            return True
    return False


def _ramp_should_suppress(t, v, intervals):
    """Return True if point (t, v) is inside a ramp window AND the flow value
    is within the normal ramp range (i.e. does not overshoot the SP target).

    Suppression rule:
      - t must be inside a ramp interval
      - If sp_target is unknown: suppress conservatively (old blanket behaviour)
      - If sp_target is known:
          suppress  when  v  <=  sp_target + tol      (tracking SP, normal ramp)
          keep scoring when  v  >  sp_target + tol      (overshoots SP, real anomaly)
        where tol = MFC_OVERSHOOT_FRAC * max(abs(sp_target), FLOW_MIN_RANGE_SCCM)
        Using FLOW_MIN_RANGE_SCCM as the floor for tol ensures that even a
        ramp-to-zero (sp_target=0) requires a meaningful overshoot before flagging.
    """
    for entry in intervals:
        t_start, t_end, sp_col, sp_target = entry
        if not (t_start <= t <= t_end):
            continue
        if sp_target is None:
            return True   # no SP info -- suppress conservatively
        tol   = max(MFC_OVERSHOOT_FRAC * abs(sp_target), FLOW_MIN_RANGE_SCCM)
        upper = sp_target + tol
        if v <= upper:
            return True   # within normal ramp range -- suppress
        return False      # significant overshoot above SP -- keep as anomaly
    return False



# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
# 4. ANOMALY DETECTION  (proven logic from old script)
# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ

def _bl_lookup(step_index, sn, se, BIN=0.5):
    """Look up (median, mad) for a given step number and step elapsed time.
    Only searches within the same step -- never crosses step boundaries.
    Returns (None, None) if the step has no baseline data.
    """
    entry = step_index.get(sn)
    if entry is None:
        return None, None
    bins = entry["bins"]
    if not bins:
        return None, None
    b = round(se / BIN) * BIN
    idx = bisect.bisect_left(bins, b)
    idx = min(idx, len(bins) - 1)
    return entry["medians"][idx], entry["mads"][idx]


def compute_run_score(rrows, baseline, sensor, step_filter=None, ramp_intervals=None):
    """exceedance-fraction x mean-excess-z for one sensor.
    Uses per-step baseline lookup -- never crosses step boundaries.
    baseline is the step_index dict returned by build_baseline.
    ramp_intervals: optional list of (t_start, t_end) to skip (MFC ramp windows).
    """
    step_index = baseline[0] if isinstance(baseline, (tuple, list)) else baseline
    if not step_index:
        return 0.0
    z_vals = []
    for r in rrows:
        v  = _f(r.get(sensor))
        se = r.get("_step_elapsed_s")
        sn = r.get("_step_number", "")
        if v is None or se is None or sn == "":
            continue
        if step_filter is not None and sn not in step_filter:
            continue
        # Skip MFC ramp windows unless the flow overshoots the SP target
        t_elapsed = r.get("_elapsed_s", 0) or 0
        if ramp_intervals and _ramp_should_suppress(t_elapsed, v, ramp_intervals):
            continue
        med, mad = _bl_lookup(step_index, sn, se)
        if med is None:
            continue
        z_vals.append(abs(v - med) / max(mad, 1e-9))
    if not z_vals:
        return 0.0
    exceed = [z for z in z_vals if z > SENSOR_DIVERGE_SIGMA]
    frac   = len(exceed) / len(z_vals)
    return frac * _mean(exceed) if exceed else 0.0
def flag_anomalous(run_composite_scores):
    """
    Adaptive threshold: median + OUTLIER_RUN_ZSCORE * max(MAD, floor).

    The old guard "if MAD/median < 20%, flag none" was designed to avoid
    false positives when ALL runs are similar -- but it fires incorrectly
    when the majority of runs score near 0 (good runs) and a minority score
    high (bad runs), because the MEDIAN of the mix is pulled toward 0 while
    the MAD is small relative to that near-zero median.

    Fix: instead of comparing MAD/median, compare the score RANGE to the
    median absolute deviation -- if any run has a score that is more than
    OUTLIER_RUN_ZSCORE * MAD above the median (with a sensible floor),
    flag it.  The old "suppress-all" branch is removed; the floor on MAD
    ensures we never flag due to pure numerical noise.

    KNOWN_ANOMALOUS runs are always flagged regardless of score.
    Returns {run_id: bool}.
    """
    scores = list(run_composite_scores.values())
    if len(scores) < 2:
        result = {rid: False for rid in run_composite_scores}
        for rid in run_composite_scores:
            if rid in KNOWN_ANOMALOUS:
                result[rid] = True
        return result

    med = _median(scores)
    mad = _median([abs(s - med) for s in scores]) * 1.4826
    threshold = med + OUTLIER_RUN_ZSCORE * max(mad, med * 0.10, 0.05)
    result = {rid: s > threshold for rid, s in run_composite_scores.items()}

    # Run-continuity check: only meaningful with enough runs (â¥10).
    # A genuine persistent fault affects consecutive runs; an isolated flagged
    # run with no adjacent flagged neighbour is likely a transient blip.
    # Demote it unless its score is dramatically above threshold (> 3Ã).
    # Skip this check for small datasets where a single anomalous run is normal.
    if len(scores) >= 10:
        try:
            run_order = sorted(run_composite_scores.keys(),
                               key=lambda x: int(x) if str(x).lstrip('-').isdigit() else x)
        except TypeError:
            run_order = sorted(run_composite_scores.keys(), key=str)
        rid_idx = {rid: i for i, rid in enumerate(run_order)}
        demoted = set()
        for rid, flagged in list(result.items()):
            if not flagged or rid in KNOWN_ANOMALOUS:
                continue
            sc = run_composite_scores[rid]
            if sc > threshold * 3:
                continue   # clearly in fault territory -- keep regardless
            i = rid_idx[rid]
            prev_ok = (i > 0 and result.get(run_order[i - 1], False))
            next_ok = (i < len(run_order) - 1 and result.get(run_order[i + 1], False))
            if not prev_ok and not next_ok:
                demoted.add(rid)
        for rid in demoted:
            result[rid] = False

    # always flag user-confirmed anomalous runs
    for rid in run_composite_scores:
        if rid in KNOWN_ANOMALOUS:
            result[rid] = True
    return result


# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
# 5. LEAD / LAG DETECTION  (ported exactly from old script)
# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ

def find_divergence_time(rrows, baseline, sensor, ramp_intervals=None):
    """
    First recipe-elapsed time (seconds) at which the sensor sustains
    MIN_CONSEC_DIVERGE consecutive points above SENSOR_DIVERGE_SIGMA,
    compared against the step-relative baseline.

    Uses a 2-second rolling median of z-scores to smooth over transient
    spikes that occur in ALL runs (e.g. TV briefly hitting 90 during
    initial pump-down before settling).  Only sustained deviations that
    persist beyond the normal transient are flagged.

    ramp_intervals: optional list of (t_start, t_end) to skip (MFC ramp windows).

    Returns recipe_elapsed_s or None.
    """
    step_index = baseline[0] if isinstance(baseline, (tuple, list)) else baseline
    if not step_index:
        return None

    t_vals, z_vals = [], []
    for r in rrows:
        v  = _f(r.get(sensor))
        se = r.get("_step_elapsed_s")
        sn = r.get("_step_number", "")
        if v is None or se is None or sn == "":
            continue
        # Skip MFC ramp windows unless the flow overshoots the SP target
        t_elapsed = r.get("_elapsed_s", 0) or 0
        if ramp_intervals and _ramp_should_suppress(t_elapsed, v, ramp_intervals):
            continue
        med, mad = _bl_lookup(step_index, sn, se)
        if med is None:
            continue
        z = abs(v - med) / max(mad, 1e-9)
        t_vals.append(r["_elapsed_s"])
        z_vals.append(z)

    if not t_vals:
        return None

    # Smooth: replace each z with the median of a +/-2s window
    SMOOTH_WIN = 2.0
    smoothed = []
    for i, t in enumerate(t_vals):
        window = [z_vals[j] for j in range(len(t_vals))
                  if abs(t_vals[j] - t) <= SMOOTH_WIN]
        smoothed.append(_median(window))

    consec = 0
    for i, z in enumerate(smoothed):
        if z > SENSOR_DIVERGE_SIGMA:
            consec += 1
            if consec >= MIN_CONSEC_DIVERGE:
                return t_vals[i - MIN_CONSEC_DIVERGE + 1]
        else:
            consec = 0
    return None


def detect_lead_lag(rrows, baselines, sensor_list, mfc_ramp_ivs=None):
    """
    Returns (lead_sensor, [(sensor, diverge_t), ...], [non_diverged...]).
    Sensors ordered by first-divergence time.

    Only considers divergence that occurs AFTER the recipe is actively running
    (elapsed_s > 1.0) to avoid false positives from different initial states
    between recipe variants.  The 2-second rolling-median smoothing in
    find_divergence_time additionally suppresses transient spikes.

    mfc_ramp_ivs: optional dict {sensor: [(t0,t1),...]} of MFC ramp intervals
    to skip during divergence detection.
    """
    # Only analyse rows from the step that triggered alarms -- determined by
    # finding the alarm step number in the run's rows.  Fall back to the
    # last step in the run if no alarm step is identified.
    # Additionally, skip the first second to avoid initial-state noise.
    step_nums = [r.get("_step_number", "") for r in rrows]
    # Find the highest step number that has sensor data (last meaningful step)
    max_step = None
    for sn in step_nums:
        try:
            v = int(sn)
            if v > 0 and (max_step is None or v > max_step):
                max_step = v
        except (ValueError, TypeError):
            pass

    if max_step is not None:
        # Only analyse the final step of the run (where the fault occurred).
        # This prevents structural differences in earlier steps (different
        # setpoints, different sub-steps) from masking the true fault sensor.
        active_rrows = [r for r in rrows
                        if _safe_int(r.get("_step_number", "")) == max_step]
        # If that's empty (e.g. step data missing), fall back to all rows > 1s
        if len(active_rrows) < 10:
            active_rrows = [r for r in rrows if r.get("_elapsed_s", 0) > 1.0]
    else:
        active_rrows = [r for r in rrows if r.get("_elapsed_s", 0) > 1.0]

    div = {}
    for s in sensor_list:
        bl = baselines.get(s)
        if bl and bl[0]:
            ivs = (mfc_ramp_ivs or {}).get(s)
            t = find_divergence_time(active_rrows, bl, s, ramp_intervals=ivs)
            if t is not None:
                div[s] = t
    if not div:
        return None, [], sensor_list[:]

    # Sort by divergence time; break ties by mean excess-z (higher z-score
    # deviation wins -- the sensor that is MORE deviated at the same onset time
    # is more likely the true cause, not a downstream effect).
    def _mean_excess_z(s):
        bl = baselines.get(s)
        if not bl or not bl[0]:
            return 0.0
        step_index = bl[0] if isinstance(bl, (tuple, list)) else bl
        zs = []
        for r in active_rrows:
            v  = _f(r.get(s))
            se = r.get("_step_elapsed_s")
            sn = r.get("_step_number", "")
            if v is None or se is None or sn == "":
                continue
            med, mad = _bl_lookup(step_index, sn, se)
            if med is None:
                continue
            zs.append(abs(v - med) / max(mad, 1e-9))
        exceed = [z for z in zs if z > SENSOR_DIVERGE_SIGMA]
        return sum(exceed) / len(exceed) if exceed else 0.0

    chain   = sorted(div.items(), key=lambda x: (x[1], -_mean_excess_z(x[0])))
    lead    = chain[0][0]
    non_div = [s for s in sensor_list if s not in div]
    return lead, chain, non_div



def _segregate_outlier_sensors(per_sensor_sc):
    """
    If a single sensor dominates scores (max score > SCORE_OUTLIER_RATIO x median
    of other sensors), return the set of such sensors so they can be scored
    independently to avoid distorting the composite threshold.

    Returns a set of sensor names to be segregated (may be empty).
    """
    # per_sensor_sc: {run_id: {sensor: score}}
    sensor_max = defaultdict(float)
    for rid, s_dict in per_sensor_sc.items():
        for s, sc in s_dict.items():
            sensor_max[s] = max(sensor_max[s], sc)

    if len(sensor_max) < 2:
        return set()

    segregated = set()
    for s, max_sc in sensor_max.items():
        # Compute median of OTHER sensors max scores
        others = sorted(sc for ss, sc in sensor_max.items() if ss != s)
        if not others:
            continue
        others_median = others[len(others) // 2]
        if others_median < 1e-6:
            continue
        if max_sc > SCORE_OUTLIER_RATIO * others_median and max_sc > 50.0:
            segregated.add(s)

    return segregated


def classify_alarm_relevance(alarm, lead_sensor, chain_sensors):
    """
    Determine whether a fault alarm is RELEVANT or IRRELEVANT to the
    lead/lag sensor cascade.

    Logic:
    1. Extract subsystem tokens from EventSource (e.g. "AT/CHA/PressCtrl"
       â ["pressctrl", "cha"]) and from EventName tokens.
    2. Extract keyword tokens from the lead and lag sensor names
       (e.g. "Chamber_Pressure" â ["chamber", "pressure"]).
    3. If ANY alarm token overlaps with ANY sensor token â RELEVANT.
    4. If the alarm is a recipe-level consequence event (e.g. CLEAN_FAILED,
       Recipe_Ended, Chamber_Changed_State) â CONSEQUENCE (not a root cause,
       but directly caused by the sensor fault).
    5. If no overlap â IRRELEVANT.

    Returns a dict:
        label    : "Relevant" | "Consequence" | "Irrelevant"
        reason   : one-line explanation string
        color    : hex color for display
    """
    event_name   = alarm.get("event_name",   "").upper()
    event_source = alarm.get("event_source", "")
    description  = alarm.get("description",  "").lower()

    # ââ Extract alarm tokens from source path and event name ââ
    # "AT/CHA/PressCtrl.PRESS_CTRL_TIMEOUT" â ["pressctrl", "press", "ctrl", "timeout", "cha"]
    def _tokens(s):
        import re
        s = s.lower()
        # split on common delimiters: / . _ space -
        parts = re.split(r'[/\._\s\-]+', s)
        return {p for p in parts if len(p) >= 3}

    alarm_tokens = _tokens(event_source) | _tokens(alarm.get("event_name", ""))

    # ââ Extract sensor keyword tokens ââ
    sensor_tokens = set()
    for s in ([lead_sensor] + list(chain_sensors)):
        sensor_tokens |= _tokens(s)

    # ââ Classification ââ

    # Recipe-level consequence alarms -- not a root cause but directly caused
    # by the sensor fault. These always follow the real alarm.
    CONSEQUENCE_KW = {
        "clean_failed", "cleanfailed", "recipe_failed", "recipefailed",
        "clean_ended", "cleanended", "recipe_ended", "recipeended",
        "chamber_changed_state", "chamberfailed", "state_changed",
        "failed", "process_failed",
    }
    event_lower = alarm.get("event_name", "").lower()
    if any(kw in event_lower for kw in CONSEQUENCE_KW) and \
       alarm.get("alarm_code", "").lower() in ("fault", "error", "") :
        # But TIMEOUT and CTRL alarms that name the sensor are the root alarm
        is_root = any(kw in event_lower for kw in
                      ["timeout", "ctrl_timeout", "not_reached", "stabilize"])
        if not is_root:
            return {
                "label":  "Consequence",
                "reason": "Recipe-level failure triggered by the sensor fault above",
                "color":  "#E65100",
            }

    # Check token overlap between alarm subsystem and sensor names
    overlap = alarm_tokens & sensor_tokens
    if overlap:
        matched = sorted(overlap)[:3]
        return {
            "label":  "Relevant",
            "reason": f"Alarm subsystem matches sensor: {', '.join(matched)}",
            "color":  "#B71C1C",
        }

    # Description text match as a fallback
    for s in ([lead_sensor] + list(chain_sensors)):
        for tok in _tokens(s):
            if len(tok) >= 5 and tok in description:
                return {
                    "label":  "Relevant",
                    "reason": f"Alarm description references '{tok}' (from {s})",
                    "color":  "#B71C1C",
                }

    return {
        "label":  "Irrelevant",
        "reason": "Alarm subsystem does not match any sensor in cascade",
        "color":  "#757575",
    }


def z_trace(rrows, baseline, sensor):
    """
    Returns (t_list, z_list) -- signed z-score at each point using the
    step-relative baseline.  X-axis is recipe_elapsed_s for display.
    """
    step_index = baseline[0] if isinstance(baseline, (tuple, list)) else baseline
    if not step_index:
        return [], []
    ts, zs = [], []
    for r in rrows:
        v  = _f(r.get(sensor))
        se = r.get("_step_elapsed_s")
        sn = r.get("_step_number", "")
        if v is None or se is None or sn == "":
            continue
        med, mad = _bl_lookup(step_index, sn, se)
        if med is None:
            continue
        z = (v - med) / max(mad, 1e-9)
        ts.append(r["_elapsed_s"])
        zs.append(z)
    return ts, zs


# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
# 6. SVG HELPERS
# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5  --  STEP IMPORTANCE, HIGHFLIERS & CORRELATION
# ══════════════════════════════════════════════════════════════════════════════

def detect_aborted_runs(tool_runs):
    """
    Tool-aware abort detection.

    Operates on the already-split tool_runs dict {tool_id: {run_id: [rows]}}
    so that runs shared across tools are evaluated independently per tool
    (uses each tool's own row count, not a cross-tool merged count).

    Returns:
      aborted    : set of (tool_id, run_id) tuples  -- composite keys
      row_counts : {(tool_id, run_id): int}          -- composite keys
    """
    row_counts: dict = {}
    aborted:    set  = set()

    for tool_id, rdict in tool_runs.items():
        if not rdict:
            continue
        tool_cnts = {rid: len(rrows) for rid, rrows in rdict.items()}
        for rid, cnt in tool_cnts.items():
            row_counts[(tool_id, rid)] = cnt

        if len(tool_cnts) < 10:
            continue   # too few runs for abort detection to be meaningful

        counts_sorted = sorted(tool_cnts.values())
        median_rows   = counts_sorted[len(counts_sorted) // 2]
        abort_thresh  = median_rows * ABORT_ROW_FRAC

        tool_aborted = {(tool_id, rid)
                        for rid, cnt in tool_cnts.items()
                        if cnt < abort_thresh}
        if tool_aborted:
            print(f'  Aborted runs for tool {tool_id} '
                  f'(rows < {abort_thresh:.0f}, median={median_rows}): '
                  f'{sorted(rid for _, rid in tool_aborted)}')
        aborted |= tool_aborted

    return aborted, row_counts


# ââ per-run process signals (arc events, alarms, step completeness) ââââââââââ

def build_run_process_signals(rows, run_col, step_col, sensor_cols):
    """
    Collect per-run signals used for highflier classification.

    Returns run_signals : { run_id: {
        'mean_dq'       : mean DATA_QUALITY value across all rows for this run
                          (None if the column is absent or all values are blank),
        'step_row_count': { step_num: row_count },
    } }
    Also returns:
        median_alarm_rows  : kept for API compatibility (always 0 now)
        median_step_rows   : { step_num: median_row_count }
    """
    dq_col = next((h for h in (rows[0].keys() if rows else [])
                   if _norm_hdr(h) in ('data_quality', 'dataquality', 'dq', 'quality')),
                  None)

    run_signals = defaultdict(lambda: {
        'mean_dq': None, 'dq_vals': [], 'step_row_count': defaultdict(int),
    })
    for row in rows:
        rid = str(row.get(run_col, '') or '').strip()
        sn  = str(row.get(step_col, '') or '').strip()
        if not rid: continue
        rs = run_signals[rid]
        if sn: rs['step_row_count'][sn] += 1
        if dq_col:
            v = _f(row.get(dq_col))
            if v is not None:
                rs['dq_vals'].append(v)

    # finalise mean_dq per run
    for rid, rs in run_signals.items():
        if rs['dq_vals']:
            rs['mean_dq'] = sum(rs['dq_vals']) / len(rs['dq_vals'])
        del rs['dq_vals']   # no longer needed

    # median step row count per step
    step_row_lists = defaultdict(list)
    for rs in run_signals.values():
        for sn, cnt in rs['step_row_count'].items():
            step_row_lists[sn].append(cnt)
    median_step_rows = {}
    for sn, cnts in step_row_lists.items():
        cnts.sort()
        median_step_rows[sn] = cnts[len(cnts) // 2]

    return dict(run_signals), 0, median_step_rows


# ââ highflier detection âââââââââââââââââââââââââââââââââââââââââââââââââââââââ

def detect_highfliers(run_means_by_step, run_order,
                      run_signals, median_alarm_rows, median_step_rows):
    """
    For each top step x sensor: flag run-means beyond median Â± 3ÃIQR.

    Classification -- two-path logic:

    PATH A  DATA_QUALITY varies across runs (IQR > 0 across per-run DQ values):
      Compute the median per-run DATA_QUALITY.  A run whose DQ is < 50% of
      that median is classified as ARTIFACT regardless of anything else.
      Runs with DQ >= 50% of median are classified as REAL.
    PATH B  DATA_QUALITY is flat (all runs identical, IQR == 0) -- the column
      carries no discriminating information.  Fall back to data-derived signals:
        ARTIFACT if run is an isolated single-sensor spike AND step rows are
                 complete (no data truncation) -- the value is a point anomaly
                 with no supporting evidence of a real fault.
        ARTIFACT if step row count < STEP_ROWS_MIN_FRAC Ã median -- the sensor
                 mean for this run is computed from an incomplete data window.
        REAL     if the run is part of a consecutive cluster of flagged runs
                 (same sensor, adjacent runs also flagged) -- sustained deviation.
        REAL     if â¥ N_SENSORS_FOR_REAL sensors are simultaneously flagged in
                 this step for this run -- multi-sensor corroboration.
        ARTIFACT otherwise (isolated single-sensor spike, data appears complete).

    In both paths the DATA_QUALITY value is noted in the reason string so the
    analyst can see what it was, even when it is not the primary classifier.

    Returns:
        highfliers : { step: { sensor: { run_id: {
            value, direction, class, reason, dq_value } } } }
    """
    # Compute per-run mean DATA_QUALITY from run_signals
    run_dq = {rid: rs.get('mean_dq', None)
              for rid, rs in run_signals.items()
              if rs.get('mean_dq') is not None}

    dq_vals = [v for v in run_dq.values() if v is not None]
    dq_iqr  = _iqr(dq_vals) if len(dq_vals) >= 4 else 0.0
    dq_med  = _median(dq_vals) if dq_vals else None
    dq_varies = dq_iqr > 0.0   # True â PATH A; False â PATH B
    dq_artifact_threshold = (dq_med * 0.50) if (dq_med and dq_varies) else None

    highfliers = {}

    for sn, rm_step in run_means_by_step.items():
        highfliers[sn] = {}
        med_step_rows = median_step_rows.get(sn, 0)

        for s, rm_sensor in rm_step.items():
            vals = list(rm_sensor.values())
            if len(vals) < 6:
                continue
            med = _median(vals)
            iqr = _iqr(vals)
            if iqr < 1e-9:
                continue
            lo = med - HIGHFLIER_IQR_MULT * iqr
            hi = med + HIGHFLIER_IQR_MULT * iqr

            flagged = {ckey: v for ckey, v in rm_sensor.items() if v < lo or v > hi}
            if not flagged:
                continue

            highfliers[sn][s] = {}
            for ckey, v in flagged.items():
                # ckey may be "tool::rid" -- extract bare rid for run_signals lookup
                bare_rid  = ckey.split("::", 1)[1] if "::" in ckey else ckey
                rs        = run_signals.get(bare_rid, {})
                idx       = run_order.index(ckey) if ckey in run_order else -1
                prev_flag = idx > 0 and run_order[idx - 1] in flagged
                next_flag = (idx >= 0 and idx < len(run_order) - 1
                             and run_order[idx + 1] in flagged)
                in_cluster    = prev_flag or next_flag
                step_rows     = rs.get('step_row_count', {}).get(sn, 0)
                step_truncated = (med_step_rows > 0
                                  and step_rows < med_step_rows * STEP_ROWS_MIN_FRAC)
                dq_val = run_dq.get(bare_rid)
                rid = ckey   # keep composite key for dict storage
                highfliers[sn][s][ckey] = dict(
                    value          = v,
                    direction      = 'HIGH' if v > hi else 'LOW',
                    dq_value       = dq_val,
                    in_cluster     = in_cluster,
                    step_truncated = step_truncated,
                    step_rows      = step_rows,
                )

        # count simultaneous sensor flags per run across all sensors this step
        run_flag_count = defaultdict(int)
        for flags in highfliers[sn].values():
            for rid in flags:
                run_flag_count[rid] += 1

        # classify each flagged (step, sensor, run) entry
        for s, flags in highfliers[sn].items():
            for rid, info in flags.items():
                dq_val     = info['dq_value']
                cluster    = info['in_cluster']
                n_sensors  = run_flag_count[rid]
                truncated  = info['step_truncated']
                step_rows  = info['step_rows']
                dq_str     = f'{dq_val:.1f}' if dq_val is not None else 'n/a'

                if dq_varies and dq_artifact_threshold is not None:
                    # ââ PATH A: DATA_QUALITY is informative âââââââââââââââ
                    if dq_val is not None and dq_val < dq_artifact_threshold:
                        cls    = 'artifact'
                        reason = (f'DATA_QUALITY = {dq_str} < 50% of median '
                                  f'({dq_med:.1f}) -- low data quality for this run')
                    else:
                        cls    = 'real'
                        reason = (f'DATA_QUALITY = {dq_str} >= 50% of median '
                                  f'({dq_med:.1f}) -- data quality acceptable; '
                                  f'deviation is a real process issue')
                else:
                    # ââ PATH B: DATA_QUALITY flat -- use data-derived signals
                    if truncated:
                        cls    = 'artifact'
                        reason = (f'Step has {step_rows} rows '
                                  f'(< {STEP_ROWS_MIN_FRAC*100:.0f}% of median '
                                  f'{med_step_rows}) -- truncated data window; '
                                  f'DQ={dq_str} (flat across all runs)')
                    elif cluster:
                        cls    = 'real'
                        reason = (f'Part of a consecutive cluster of flagged runs -- '
                                  f'sustained deviation; '
                                  f'DQ={dq_str} (flat, not discriminating)')
                    elif n_sensors >= N_SENSORS_FOR_REAL:
                        cls    = 'real'
                        reason = (f'{n_sensors} sensors simultaneously elevated -- '
                                  f'multi-sensor corroboration; '
                                  f'DQ={dq_str} (flat, not discriminating)')
                    else:
                        cls    = 'artifact'
                        reason = (f'Isolated single-sensor spike, step data complete '
                                  f'({step_rows} rows); '
                                  f'DQ={dq_str} (flat, not discriminating)')

                info['class']  = cls
                info['reason'] = reason

    return highfliers


# ââ correlation matrix ââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
def compute_correlations(rm_step, exclude_runs):
    """
    Compute Pearson correlation matrix of per-run means across sensors.
    Excludes aborted runs and runs classified as data artifacts.

    rm_step      : { sensor: { run_id: mean_value } }
    exclude_runs : set of run_ids to exclude

    Returns (sensors, corr_matrix) where corr_matrix[i][j] is the Pearson r
    between sensor i and sensor j.  Returns ([], []) if < 3 usable runs.

    Sensors that have no data at all for the non-excluded runs are dropped
    before computing the intersection so that tool-specific missing sensors
    (e.g. a column that is always null for one tool) do not zero out the
    common-run set for that tool.
    """
    # Collect sensors that have at least one non-excluded run
    all_sensors = sorted(rm_step.keys())
    if len(all_sensors) < 2:
        return [], []

    # First pass: find which runs are available after exclusion, per sensor
    sensor_runs = {s: set(rm_step[s].keys()) - exclude_runs for s in all_sensors}

    # Drop sensors that have zero non-excluded runs for this tool subset --
    # these are columns missing entirely for the target tool and would make
    # the intersection empty, hiding all correlations for that tool.
    sensors = [s for s in all_sensors if len(sensor_runs[s]) >= 3]
    if len(sensors) < 2:
        return [], []

    # find runs present in ALL remaining sensors and not excluded
    common = sensor_runs[sensors[0]]
    for s in sensors[1:]:
        common &= sensor_runs[s]
    runs = sorted(common)
    if len(runs) < 3:
        return [], []

    # build matrix: values[i] = list of per-run means for sensor i
    values = [[rm_step[s][rid] for rid in runs] for s in sensors]
    n      = len(runs)

    def _pearson(x, y):
        mx = sum(x)/n; my = sum(y)/n
        num  = sum((x[k]-mx)*(y[k]-my) for k in range(n))
        d1   = math.sqrt(sum((x[k]-mx)**2 for k in range(n)))
        d2   = math.sqrt(sum((y[k]-my)**2 for k in range(n)))
        return num / max(d1*d2, 1e-12)

    ns = len(sensors)
    corr = [[_pearson(values[i], values[j]) for j in range(ns)] for i in range(ns)]
    return sensors, corr

def score_steps(rows, run_col, tool_col, step_col, sensor_cols, aborted_runs,
                noisy_sensor_steps=None):
    """Score steps by cross-run CV, excluding aborted runs.

    Uses composite key  tool + '::' + rid  so that the same run ID appearing
    on multiple tools is tracked independently.  run_order and run_tool use
    these composite keys; callers receive them and pass them on to
    collect_run_means / svg_trend_chart unchanged.

    noisy_sensor_steps: set of (sensor, step_key) pairs to exclude from
      scoring.  Populated by filter_sensor_cols from the per-step noise analysis.
    """
    noisy_sensor_steps = noisy_sensor_steps or set()
    accum           = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    step_row_counts = defaultdict(int)
    step_names      = {}
    run_order       = []; seen_runs = set()
    run_tool        = {}

    step_name_col = None
    if rows:
        for h in rows[0].keys():
            if _norm_hdr(h) in ('step_name','stepname','step_description'):
                step_name_col = h; break

    for row in rows:
        rid  = str(row.get(run_col,  '') or '').strip()
        sn   = str(row.get(step_col, '') or '').strip()
        tool = str(row.get(tool_col, '') or 'ALL').strip() if tool_col else 'ALL'
        if not rid or not sn or sn in SKIP_STEPS:
            continue
        # Composite key: tool::rid so same run ID on different tools is distinct
        ckey = f"{tool}::{rid}"
        if ckey not in seen_runs:
            seen_runs.add(ckey); run_order.append(ckey)
            run_tool[ckey] = tool
        if rid in aborted_runs:
            continue    # exclude aborted runs from all step scoring
        step_row_counts[sn] += 1
        if step_name_col and sn not in step_names:
            nm = str(row.get(step_name_col, '') or '').strip()
            if nm: step_names[sn] = nm
        for s in sensor_cols:
            v = _f(row.get(s))
            if v is not None:
                accum[sn][s][ckey].append(v)

    # Global min/max per sensor across all steps (active-step filter)
    global_sensor_min = {}; global_sensor_max = {}
    for sn, sens_dict in accum.items():
        for s, run_dict in sens_dict.items():
            for vals in run_dict.values():
                if vals:
                    lo = min(vals); hi = max(vals)
                    global_sensor_min[s] = min(global_sensor_min.get(s, lo), lo)
                    global_sensor_max[s] = max(global_sensor_max.get(s, hi), hi)

    step_info = {}
    for sn, sens_dict in accum.items():
        if step_row_counts[sn] < MIN_STEP_ROWS:
            continue
        sensor_cvs = []
        for s, run_dict in sens_dict.items():
            # Skip (sensor, step) pairs identified as noisy in filter_sensor_cols
            if (s, sn) in noisy_sensor_steps:
                continue
            run_means = [sum(v)/len(v) for v in run_dict.values() if v]
            if len(run_means) < MIN_RUNS_FOR_SCORE:
                continue
            mn  = sum(run_means)/len(run_means)
            if abs(mn) < 1e-9:
                continue
            std = math.sqrt(sum((v-mn)**2 for v in run_means)/len(run_means))
            cv  = std / abs(mn)
            rng = max(run_means) - min(run_means)

            g_lo  = global_sensor_min.get(s, mn)
            g_hi  = global_sensor_max.get(s, mn)
            g_rng = max(g_hi - g_lo, 1e-12)

            # Use the step-level range of run-means (rng) as the reference for
            # the minimum-std check.  Using the global sensor range (g_rng) is
            # too strict for sensors whose step-level variation is meaningful but
            # small relative to their full operating range across all steps
            # (e.g. Pedestal heater varies 0.5% at Dep but spans 5–60% globally).
            step_rng = max(rng, 1e-12)
            if std < step_rng * MIN_STD_FRAC_OF_RANGE:
                continue
            # Skip sensors whose run-means band is too tight to be visually
            # meaningful -- the full spread of run means must be >= 3% of the
            # step mean.  This excludes mechanically-locked sensors (e.g.
            # Pedestal lift position always within 1.86 mm of 84 mm = 2.2%).
            if abs(mn) > 1e-9 and rng < MIN_RNG_FRAC_OF_MEAN * abs(mn):
                continue
            # Skip sensors whose mean is at/near the inactive minimum.
            # For sensors that never reach zero (always-negative), skip this
            # check as they are always "active" from a process perspective.
            active_threshold = g_lo + g_rng * MIN_ACTIVE_FRAC
            if g_lo >= 0 and mn < active_threshold:
                continue
            elif g_lo < 0 and g_hi < 0:
                pass  # always-negative sensor: always active, skip the check
            elif mn < active_threshold:
                continue
            # Late-shift detection: boost CV if sensor shows a significant
            # trend shift between the first and last third of runs.
            run_means_sorted = list(run_means)
            n_early = max(1, len(run_means_sorted) // 3)
            n_late  = max(1, len(run_means_sorted) // 3)
            early_mean = sum(run_means_sorted[:n_early]) / n_early
            late_mean  = sum(run_means_sorted[-n_late:]) / n_late
            shift_magnitude = abs(late_mean - early_mean)
            shift_score = shift_magnitude / std if std > 0 else 0.0
            # Boost effective CV if there is a late-run shift
            effective_cv = max(cv, cv * (1 + shift_score * 0.5))
            if effective_cv >= MIN_CV_FOR_CHART:
                sensor_cvs.append((s, effective_cv, mn, std, rng))

        if not sensor_cvs:
            continue
        sensor_cvs.sort(key=lambda x: -x[1])
        mean_cv = sum(x[1] for x in sensor_cvs) / len(sensor_cvs)
        n_runs  = max(len(rd) for rd in sens_dict.values())
        step_info[sn] = dict(
            name       = step_names.get(sn, f'Step {sn}'),
            rows       = step_row_counts[sn],
            n_runs     = n_runs,
            mean_cv    = mean_cv,
            sensor_cvs = sensor_cvs,
        )

    return step_info, run_order, run_tool, step_names


# ââ per-run means for trend charts âââââââââââââââââââââââââââââââââââââââââââ
def collect_run_means(rows, run_col, step_col, sensors, target_steps, aborted_runs,
                      tool_col=None):
    """Returns { step: { sensor: { tool::rid: mean_value } } } -- aborted excluded.

    Uses the same tool::rid composite key as score_steps so that shared run IDs
    across multiple tools are tracked independently in trend charts.
    """
    accum = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for row in rows:
        rid  = str(row.get(run_col,  '') or '').strip()
        sn   = str(row.get(step_col, '') or '').strip()
        tool = str(row.get(tool_col, '') or 'ALL').strip() if tool_col else 'ALL'
        if not rid or rid in aborted_runs or sn not in target_steps:
            continue
        ckey = f"{tool}::{rid}"
        for s in sensors:
            v = _f(row.get(s))
            if v is not None:
                accum[sn][s][ckey].append(v)
    result = {}
    for sn, sd in accum.items():
        result[sn] = {}
        for s, rd in sd.items():
            result[sn][s] = {ckey: sum(v)/len(v) for ckey, v in rd.items() if v}
    return result


# ââ SVG trend chart âââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ

TOOL_COLORS = [
    '#1565C0','#2E7D32','#6A1B9A','#C62828','#E65100',
    '#00838F','#4E342E','#37474F','#AD1457','#558B2F',
]
ABORT_COLOR      = '#E65100'   # orange triangle for aborted runs
HIGHFLIER_REAL   = '#B71C1C'   # dark-red diamond -- real issue
HIGHFLIER_ART    = '#FBC02D'   # amber diamond -- data artifact

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6  --  SVG CHART FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def _esc(s):
    return (str(s).replace("&","&amp;").replace("<","&lt;")
                  .replace(">","&gt;").replace('"',"&quot;"))

def _fmt_t(t):
    if t < 60:
        return f"{t:.0f}s"
    m = int(t // 60)
    s = t % 60
    return f"{m}m {s:.0f}s" if s else f"{m}m"

def _fmt_v(v):
    """Format a sensor value as plain decimal -- never use exponential notation."""
    if v == 0:
        return "0"
    abs_v = abs(v)
    if abs_v >= 1000:
        return f"{v:.0f}"
    if abs_v >= 100:
        return f"{v:.1f}"
    if abs_v >= 10:
        return f"{v:.2f}"
    if abs_v >= 1:
        return f"{v:.3f}"
    if abs_v >= 0.01:
        return f"{v:.4f}"
    # very small: use enough decimal places to show significant digits
    return f"{v:.6f}".rstrip('0').rstrip('.')

def _nice_interval(t_range, max_ticks=12):
    raw = t_range / max_ticks
    for step in [1,2,5,10,15,20,30,60,120,300,600]:
        if step >= raw:
            return step
    return int(raw) + 1


# ââ A. Raw sensor trace overlay âââââââââââââââââââââââââââââââââââââââââââââââ
def svg_raw_trace(sensor, tool_id, run_dict, step_info,
                  looped_per_run, run_colors, anomalous_runs,
                  alarm_info=None, aborted_runs=None):
    """
    All runs overlaid for one sensor with:
    - Shaded envelope (min/max of normal runs) as a green band
    - Normal runs as thin colored solid lines
    - Anomalous run as a thick bright-red line (clearly distinct)
    - Step bands with legible labels in a dedicated label row below X-axis
    - Legend outside the plot area on the right
    - Alarm markers as bold vertical lines with text annotation boxes
    """
    traces = {}
    for rid, rrows in run_dict.items():
        pts_t, pts_v = [], []
        for r in rrows:
            v = _f(r.get(sensor))
            if v is not None:
                pts_t.append(r["_elapsed_s"])
                pts_v.append(v)
        if pts_t:
            traces[rid] = (pts_t, pts_v)
    if not traces:
        return ""

    # reference run = most steps
    ref    = max(run_dict, key=lambda r: len(step_info.get(r, [])))
    all_ts = [t for ts, _ in traces.values() for t in ts]
    t0, tmax = min(all_ts), max(all_ts)

    # deduplicate step boundaries (same step# + name = same step)
    seen, clean = set(), []
    for b in step_info.get(ref, []):
        key = (b["step_number"], b["step_name"])
        if key in seen or b["elapsed_s"] > tmax + 5:
            continue
        seen.add(key)
        clean.append(b)
    n_steps = len(clean)

    # layout -- legend on the right, step labels below in a dedicated strip
    LEG_W    = 220   # right legend panel width
    LABEL_H  = 90    # bottom strip height for step labels
    mg = {"top": 65, "right": LEG_W + 20, "bottom": LABEL_H, "left": 85}
    W  = max(1000, min(MAX_CHART_WIDTH, n_steps * MIN_PX_PER_STEP + mg["left"] + mg["right"]))
    H  = CHART_HEIGHT
    pw = W - mg["left"] - mg["right"]
    ph = H - mg["top"]  - mg["bottom"]

    # Setup aborted run set and normal_ids (aborted excluded from envelope)
    _aborted_set = aborted_runs or set()
    # value range -- from normal runs only (so anomalous run sticks out visually)
    normal_ids = [r for r in traces if r not in anomalous_runs and r not in _aborted_set]
    if normal_ids:
        norm_vs = [v for rid in normal_ids for v in traces[rid][1]]
    else:
        norm_vs = [v for _, vs in traces.values() for v in vs]
    all_vs   = [v for _, vs in traces.values() for v in vs]
    vmin_data, vmax_data = min(all_vs), max(all_vs)
    pad  = (vmax_data - vmin_data) * 0.12 if vmax_data != vmin_data else 1.0
    vmin = vmin_data - pad
    vmax = vmax_data + pad
    tr   = max(tmax - t0, 1e-9)
    vr   = max(vmax - vmin, 1e-9)

    def tx(t): return mg["left"] + (t - t0) / tr * pw
    def ty(v): return mg["top"]  + (1 - (v - vmin) / vr) * ph

    svg = []
    svg.append(f'<svg width="{W}" height="{H}" xmlns="http://www.w3.org/2000/svg" '
               f'style="background:#FAFAFA;font-family:Arial,sans-serif;display:block;margin:10px 0">')

    # ââ title ââ
    svg.append(f'<text x="{mg["left"] + pw//2}" y="22" text-anchor="middle" '
               f'font-size="15" font-weight="bold" fill="#1A237E">'
               f'Tool {_esc(tool_id)} -- {_esc(sensor)}</text>')
    anom_label = (", ".join(f"Run {r}" for r in sorted(anomalous_runs) if r in traces)
                  or "none")
    svg.append(f'<text x="{mg["left"] + pw//2}" y="40" text-anchor="middle" '
               f'font-size="11" fill="#555">'
               f'{n_steps} steps  |  {len(traces)} runs overlaid  |  '
               f'Anomalous: {_esc(anom_label)}</text>')

    # ââ plot area background ââ
    svg.append(f'<rect x="{mg["left"]}" y="{mg["top"]}" width="{pw}" height="{ph}" '
               f'fill="white" stroke="#B0BEC5" stroke-width="1.5"/>')

    # ââ horizontal grid lines ââ
    for i in range(7):
        gv = vmin + i * vr / 6
        gy = ty(gv)
        svg.append(f'<line x1="{mg["left"]}" y1="{gy:.1f}" '
                   f'x2="{mg["left"]+pw}" y2="{gy:.1f}" '
                   f'stroke="#ECEFF1" stroke-width="1"/>')

    # ââ step background bands (drawn before traces) ââ
    looped = looped_per_run.get(ref, set())
    for bi, b in enumerate(clean):
        bx = tx(b["elapsed_s"])
        ex = tx(clean[bi + 1]["elapsed_s"]) if bi + 1 < len(clean) else tx(tmax)
        sw = max(ex - bx, 1)
        is_lp = b["step_name"] in looped
        if is_lp:
            fill = "#FFF8E1" if bi % 2 == 0 else "#FFF3E0"
        else:
            fill = "white" if bi % 2 == 0 else "#F5F5F5"
        svg.append(f'<rect x="{bx:.1f}" y="{mg["top"]}" width="{sw:.1f}" '
                   f'height="{ph}" fill="{fill}" stroke="none"/>')
        # solid thin line at step boundary
        lc = "#FF8F00" if is_lp else "#90A4AE"
        lw = "1.5" if is_lp else "1"
        svg.append(f'<line x1="{bx:.1f}" y1="{mg["top"]}" '
                   f'x2="{bx:.1f}" y2="{mg["top"]+ph}" '
                   f'stroke="{lc}" stroke-width="{lw}" opacity="0.7"/>')

    # ââ envelope band: min/max of normal runs at each time bin ââ
    # Build a time-sorted list of (t, min_v, max_v) across normal runs
    if len(normal_ids) >= 2:
        BIN = 1.0   # 1-second bins for envelope
        env = defaultdict(list)
        for rid in normal_ids:
            for t, v in zip(*traces[rid]):
                env[round(t / BIN) * BIN].append(v)
        env_times  = sorted(env)
        env_lo     = [min(env[t]) for t in env_times]
        env_hi     = [max(env[t]) for t in env_times]
        # polygon: upper edge LâR then lower edge RâL
        upper = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(env_times, env_hi))
        lower = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(reversed(env_times), reversed(env_lo)))
        svg.append(f'<polygon points="{upper} {lower}" '
                   f'fill="#A5D6A7" opacity="0.35" stroke="none"/>')
        # envelope border lines
        upper_pts = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(env_times, env_hi))
        lower_pts = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(env_times, env_lo))
        svg.append(f'<polyline points="{upper_pts}" fill="none" '
                   f'stroke="#388E3C" stroke-width="1" opacity="0.5"/>')
        svg.append(f'<polyline points="{lower_pts}" fill="none" '
                   f'stroke="#388E3C" stroke-width="1" opacity="0.5"/>')

    # ââ normal run traces (thin, semi-transparent) ââ
    for rid in sorted(traces):
        if rid in anomalous_runs or rid in _aborted_set:
            continue
        ts, vs = traces[rid]
        col  = run_colors.get(rid, "#777")
        pts  = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(ts, vs))
        svg.append(f'<polyline points="{pts}" fill="none" '
                   f'stroke="{col}" stroke-width="1.8" opacity="0.65"/>')

    # aborted run traces (orange, medium weight)
    if _aborted_set:
        for rid in sorted(_aborted_set):
            if rid not in traces:
                continue
            ts, vs = traces[rid]
            pts = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(ts, vs))
            svg.append(f'<polyline points="{pts}" fill="none" '
                       f'stroke="#FF6D00" stroke-width="2" opacity="0.75"/>')

    # anomalous run traces (thick, bright red, on top)
    for rid in sorted(anomalous_runs):
        if rid not in traces:
            continue
        ts, vs = traces[rid]
        pts = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(ts, vs))
        # thick white halo first (makes it stand out on any background)
        svg.append(f'<polyline points="{pts}" fill="none" '
                   f'stroke="white" stroke-width="6" opacity="0.7"/>')
        svg.append(f'<polyline points="{pts}" fill="none" '
                   f'stroke="#D32F2F" stroke-width="3.5" opacity="1.0"/>')

    # ââ alarm markers ââ
    # ARC onset events â purple diamond + dashed line
    # Regular alarms  â red diamond + dashed line
    alarm_annotations = []   # collect for text boxes below chart
    if alarm_info:
        for rid in sorted(anomalous_runs):
            for a in (alarm_info.get(rid) or []):
                at = a["elapsed_s"]
                if at < t0 or at > tmax:
                    continue
                ax = tx(at)
                is_arc = a.get("is_arc_onset") or a.get("alarm_code") == "ARC_EVENT"
                mc = "#6A1B9A" if is_arc else "#B71C1C"
                svg.append(f'<line x1="{ax:.1f}" y1="{mg["top"]}" '
                           f'x2="{ax:.1f}" y2="{mg["top"]+ph}" '
                           f'stroke="{mc}" stroke-width="2.5" '
                           f'stroke-dasharray="8,4" opacity="0.95"/>')
                dy = mg["top"] - 6
                svg.append(f'<polygon points="{ax:.1f},{dy-8} {ax+6:.1f},{dy} '
                           f'{ax:.1f},{dy+8} {ax-6:.1f},{dy}" '
                           f'fill="{mc}"/>')
                alarm_annotations.append((ax, a))

    # ââ Y-axis ticks and label ââ
    for i in range(7):
        v  = vmin + i * vr / 6
        yp = ty(v)
        svg.append(f'<line x1="{mg["left"]-6}" y1="{yp:.1f}" '
                   f'x2="{mg["left"]}" y2="{yp:.1f}" '
                   f'stroke="#90A4AE" stroke-width="1.5"/>')
        svg.append(f'<text x="{mg["left"]-10}" y="{yp+4:.1f}" '
                   f'text-anchor="end" font-size="11" fill="#37474F">'
                   f'{_fmt_v(v)}</text>')
    midy = mg["top"] + ph / 2
    svg.append(f'<text x="16" y="{midy:.1f}" text-anchor="middle" '
               f'font-size="12" fill="#455A64" '
               f'transform="rotate(-90,16,{midy:.1f})">{_esc(sensor)}</text>')

    # ââ X-axis ticks ââ
    ti = _nice_interval(tr)
    tt = 0.0
    while tt <= tmax + 0.5:
        xp = tx(tt)
        svg.append(f'<line x1="{xp:.1f}" y1="{mg["top"]+ph}" '
                   f'x2="{xp:.1f}" y2="{mg["top"]+ph+6}" '
                   f'stroke="#90A4AE" stroke-width="1.5"/>')
        svg.append(f'<text x="{xp:.1f}" y="{mg["top"]+ph+20}" '
                   f'text-anchor="middle" font-size="11" fill="#37474F">'
                   f'{_fmt_t(tt)}</text>')
        tt += ti
    svg.append(f'<text x="{mg["left"]+pw//2}" y="{mg["top"]+ph+38}" '
               f'text-anchor="middle" font-size="12" fill="#455A64" font-weight="bold">'
               f'Recipe Elapsed Time (s)</text>')

    # ââ step labels in dedicated label strip below X ticks ââ
    label_y0 = mg["top"] + ph + 52   # top of label strip
    for bi, b in enumerate(clean):
        bx = tx(b["elapsed_s"])
        is_lp = b["step_name"] in looped
        lc = "#BF360C" if is_lp else "#37474F"
        fw = "bold" if is_lp else "normal"
        lbl = f'{b["step_number"]}, {(b["step_name"] or "?")[:16]}'
        svg.append(f'<text x="{bx+3:.1f}" y="{label_y0}" '
                   f'font-size="9" fill="{lc}" font-weight="{fw}" '
                   f'transform="rotate(55,{bx+3:.1f},{label_y0})">'
                   f'{_esc(lbl)}</text>')

    # ââ right-side legend panel ââ
    lx = mg["left"] + pw + 18
    ly = mg["top"]
    row = 0

    def leg_row(y_off, line_color, line_w, line_dash, text, text_color="#222", bold=False):
        ey = ly + y_off + 8
        d  = f' stroke-dasharray="{line_dash}"' if line_dash else ""
        svg.append(f'<line x1="{lx}" y1="{ey}" x2="{lx+28}" y2="{ey}" '
                   f'stroke="{line_color}" stroke-width="{line_w}"{d}/>')
        fw2 = "bold" if bold else "normal"
        svg.append(f'<text x="{lx+34}" y="{ey+4}" font-size="11" '
                   f'fill="{text_color}" font-weight="{fw2}">{_esc(text)}</text>')

    # envelope swatch
    svg.append(f'<rect x="{lx}" y="{ly+4}" width="28" height="12" '
               f'fill="#A5D6A7" opacity="0.5" stroke="#388E3C" stroke-width="1"/>')
    svg.append(f'<text x="{lx+34}" y="{ly+14}" font-size="11" fill="#2E7D32">'
               f'Normal envelope</text>')
    row_off = 24

    # normal runs
    for ri, rid in enumerate(sorted(r for r in traces if r not in anomalous_runs and r not in _aborted_set)):
        col = run_colors.get(rid, "#777")
        leg_row(row_off, col, "1.8", "", f"Run {rid} (normal)", col)
        row_off += 20

    # aborted runs legend
    for rid in sorted(_aborted_set & set(traces)):
        leg_row(row_off, "#FF6D00", "2", "", f"Run {rid}  ABORTED", "#FF6D00", bold=False)
        row_off += 20

    # anomalous runs
    for rid in sorted(anomalous_runs & set(traces)):
        leg_row(row_off, "#D32F2F", "3.5", "", f"Run {rid}  ANOMALOUS", "#C62828", bold=True)
        row_off += 20

    # alarm marker
    if alarm_annotations:
        ey2 = ly + row_off + 8
        svg.append(f'<line x1="{lx}" y1="{ey2}" x2="{lx+28}" y2="{ey2}" '
                   f'stroke="#B71C1C" stroke-width="2.5" stroke-dasharray="8,4"/>')
        svg.append(f'<polygon points="{lx+14},{ey2-10} {lx+20},{ey2} '
                   f'{lx+14},{ey2+10} {lx+8},{ey2}" fill="#B71C1C"/>')
        svg.append(f'<text x="{lx+34}" y="{ey2+4}" font-size="11" '
                   f'fill="#B71C1C" font-weight="bold">Fault alarm</text>')
        row_off += 28
        # alarm descriptions
        for ax_px, a in alarm_annotations:
            row_off += 2
            short = a["event_name"].split(".")[-1][:30]
            desc  = a["description"][:50].strip()
            svg.append(f'<text x="{lx}" y="{ly+row_off+10}" font-size="9" '
                       f'fill="#B71C1C" font-weight="bold">{_esc(short)}</text>')
            row_off += 14
            svg.append(f'<text x="{lx}" y="{ly+row_off+10}" font-size="9" '
                       f'fill="#555">{_esc(desc)}</text>')
            row_off += 16

    # step type legend
    row_off += 8
    svg.append(f'<rect x="{lx}" y="{ly+row_off}" width="14" height="12" '
               f'fill="#FFF3E0" stroke="#FF8F00" stroke-width="1"/>')
    svg.append(f'<text x="{lx+20}" y="{ly+row_off+10}" font-size="10" '
               f'fill="#BF360C">Looped step</text>')
    row_off += 18
    svg.append(f'<rect x="{lx}" y="{ly+row_off}" width="14" height="12" '
               f'fill="#F5F5F5" stroke="#90A4AE" stroke-width="1"/>')
    svg.append(f'<text x="{lx+20}" y="{ly+row_off+10}" font-size="10" '
               f'fill="#546E7A">Normal step</text>')

    svg.append("</svg>")
    return "\n".join(svg)


# ââ B2. Baseline envelope trace chart âââââââââââââââââââââââââââââââââââââââââ

def svg_envelope_trace(sensor, tool_id, run_dict, step_info,
                       looped_per_run, run_colors, anomalous_runs,
                       alarm_info=None, aborted_runs=None):
    """
    Baseline envelope chart for one sensor -- mirrors the reference good-vs-bad
    images:

      â¢ Green shaded band  = median Â± 3ÃMAD of normal runs (at each time bin)
      â¢ Thin green lines   = individual normal runs
      â¢ Thick red lines    = anomalous runs (white halo underneath for contrast)
      â¢ Vertical dashed lines = step boundaries
      â¢ Step labels below X-axis
      â¢ Legend on the right

    Unlike svg_raw_trace (which scales the Y-axis from normal runs only and
    uses the raw min/max envelope), this chart uses a statistics-based band
    so the envelope is tight around the typical behaviour and anomalous
    excursions stand out clearly -- exactly as in the reference images.
    """
    traces = {}
    for rid, rrows in run_dict.items():
        pts_t, pts_v = [], []
        for r in rrows:
            v = _f(r.get(sensor))
            if v is not None:
                pts_t.append(r["_elapsed_s"])
                pts_v.append(v)
        if pts_t:
            traces[rid] = (pts_t, pts_v)
    if not traces:
        return ""

    _aborted_set_env = aborted_runs or set()
    normal_ids = [r for r in traces if r not in anomalous_runs and r not in _aborted_set_env]
    anom_ids   = [r for r in traces if r in anomalous_runs]

    # reference run = most steps
    ref    = max(run_dict, key=lambda r: len(step_info.get(r, [])))
    all_ts = [t for ts, _ in traces.values() for t in ts]
    t0, tmax = min(all_ts), max(all_ts)

    # deduplicate step boundaries
    seen, clean = set(), []
    for b in step_info.get(ref, []):
        key = (b["step_number"], b["step_name"])
        if key in seen or b["elapsed_s"] > tmax + 5:
            continue
        seen.add(key)
        clean.append(b)
    n_steps = len(clean)

    # layout
    LEG_W   = 230
    LABEL_H = 90
    mg = {"top": 65, "right": LEG_W + 20, "bottom": LABEL_H, "left": 85}
    W  = max(1000, min(MAX_CHART_WIDTH, n_steps * MIN_PX_PER_STEP + mg["left"] + mg["right"]))
    H  = CHART_HEIGHT
    pw = W - mg["left"] - mg["right"]
    ph = H - mg["top"]  - mg["bottom"]

    # ââ build statistics-based envelope from normal runs ââ
    # Bin data at 1-second resolution across all normal runs
    BIN = 1.0
    env_bins = defaultdict(list)
    for rid in normal_ids:
        for t, v in zip(*traces[rid]):
            env_bins[round(t / BIN) * BIN].append(v)

    env_times, env_lo, env_hi, env_med = [], [], [], []
    for tb in sorted(env_bins):
        vs = env_bins[tb]
        if not vs:
            continue
        med = sorted(vs)[len(vs) // 2]
        mad = sorted(abs(v - med) for v in vs)[len(vs) // 2] if len(vs) > 1 else 0.0
        env_times.append(tb)
        env_med.append(med)
        env_lo.append(med - 3.0 * mad)
        env_hi.append(med + 3.0 * mad)

    # Y range: accommodate both envelope and all anomalous traces
    all_vs = [v for _, vs in traces.values() for v in vs]
    lo_data = min(all_vs) if all_vs else 0
    hi_data = max(all_vs) if all_vs else 1
    lo_env  = min(env_lo)  if env_lo  else lo_data
    hi_env  = max(env_hi)  if env_hi  else hi_data
    vmin_data = min(lo_data, lo_env)
    vmax_data = max(hi_data, hi_env)
    pad  = (vmax_data - vmin_data) * 0.12 if vmax_data != vmin_data else 1.0
    vmin = vmin_data - pad
    vmax = vmax_data + pad
    tr   = max(tmax - t0, 1e-9)
    vr   = max(vmax - vmin, 1e-9)

    def tx(t): return mg["left"] + (t - t0) / tr * pw
    def ty(v): return mg["top"]  + (1 - (v - vmin) / vr) * ph

    svg = []
    svg.append(f'<svg width="{W}" height="{H}" xmlns="http://www.w3.org/2000/svg" '
               f'style="background:#FAFAFA;font-family:Arial,sans-serif;display:block;margin:10px 0">')

    # title
    anom_label = ", ".join(f"Run {r}" for r in sorted(anom_ids)) or "none"
    svg.append(f'<text x="{mg["left"] + pw//2}" y="22" text-anchor="middle" '
               f'font-size="15" font-weight="bold" fill="#1A237E">'
               f'Tool {_esc(tool_id)} -- {_esc(sensor)}</text>')
    svg.append(f'<text x="{mg["left"] + pw//2}" y="40" text-anchor="middle" '
               f'font-size="11" fill="#555">'
               f'Baseline envelope (median Â± 3ÃMAD of {len(normal_ids)} normal runs)  |  '
               f'Anomalous: {_esc(anom_label)}</text>')

    # plot area
    svg.append(f'<rect x="{mg["left"]}" y="{mg["top"]}" width="{pw}" height="{ph}" '
               f'fill="white" stroke="#B0BEC5" stroke-width="1.5"/>')

    # horizontal grid lines
    for i in range(7):
        gv = vmin + i * vr / 6
        gy = ty(gv)
        svg.append(f'<line x1="{mg["left"]}" y1="{gy:.1f}" '
                   f'x2="{mg["left"]+pw}" y2="{gy:.1f}" '
                   f'stroke="#ECEFF1" stroke-width="1"/>')

    # step bands and boundary lines (dashed)
    looped = looped_per_run.get(ref, set())
    for bi, b in enumerate(clean):
        bx = tx(b["elapsed_s"])
        ex = tx(clean[bi + 1]["elapsed_s"]) if bi + 1 < len(clean) else tx(tmax)
        sw = max(ex - bx, 1)
        is_lp = b["step_name"] in looped
        fill = ("#FFF8E1" if bi % 2 == 0 else "#FFF3E0") if is_lp else ("white" if bi % 2 == 0 else "#F5F5F5")
        svg.append(f'<rect x="{bx:.1f}" y="{mg["top"]}" width="{sw:.1f}" '
                   f'height="{ph}" fill="{fill}" stroke="none"/>')
        lc = "#FF8F00" if is_lp else "#90A4AE"
        svg.append(f'<line x1="{bx:.1f}" y1="{mg["top"]}" '
                   f'x2="{bx:.1f}" y2="{mg["top"]+ph}" '
                   f'stroke="{lc}" stroke-width="1" stroke-dasharray="4,3" opacity="0.7"/>')

    # ââ shaded envelope band (green) ââ
    if len(env_times) >= 2:
        upper = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(env_times, env_hi))
        lower = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(reversed(env_times), reversed(env_lo)))
        svg.append(f'<polygon points="{upper} {lower}" '
                   f'fill="#A5D6A7" opacity="0.40" stroke="none"/>')
        # envelope border lines
        svg.append(f'<polyline points="{upper}" fill="none" '
                   f'stroke="#2E7D32" stroke-width="1.2" opacity="0.6"/>')
        lower_fwd = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(env_times, env_lo))
        svg.append(f'<polyline points="{lower_fwd}" fill="none" '
                   f'stroke="#2E7D32" stroke-width="1.2" opacity="0.6"/>')
        # median line
        med_pts = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(env_times, env_med))
        svg.append(f'<polyline points="{med_pts}" fill="none" '
                   f'stroke="#1B5E20" stroke-width="1.5" stroke-dasharray="6,3" opacity="0.7"/>')

    # ââ normal run traces (thin green) ââ
    for rid in sorted(normal_ids):
        ts, vs = traces[rid]
        col = run_colors.get(rid, "#388E3C")
        pts = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(ts, vs))
        svg.append(f'<polyline points="{pts}" fill="none" '
                   f'stroke="{col}" stroke-width="1.2" opacity="0.50"/>')

    # ââ anomalous run traces (thick red, white halo) ââ
    # aborted run traces (orange, medium weight)
    if _aborted_set_env:
        for rid in sorted(_aborted_set_env):
            if rid not in traces:
                continue
            ts, vs = traces[rid]
            pts = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(ts, vs))
            svg.append(f'<polyline points="{pts}" fill="none" '
                       f'stroke="#FF6D00" stroke-width="2" opacity="0.75"/>')

    # anomalous run traces (thick red, white halo)
    for rid in sorted(anom_ids):
        ts, vs = traces[rid]
        pts = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(ts, vs))
        svg.append(f'<polyline points="{pts}" fill="none" '
                   f'stroke="white" stroke-width="6" opacity="0.75"/>')
        svg.append(f'<polyline points="{pts}" fill="none" '
                   f'stroke="#D32F2F" stroke-width="3.0" opacity="1.0"/>')

    # ââ alarm markers for anomalous runs ââ
    # ARC onset â purple; regular alarms â red
    alarm_annotations = []
    for rid in sorted(anom_ids):
        for a in (alarm_info or {}).get(rid, []):
            ax_px = tx(a["elapsed_s"])
            if mg["left"] <= ax_px <= mg["left"] + pw:
                is_arc = a.get("is_arc_onset") or a.get("alarm_code") == "ARC_EVENT"
                mc = "#6A1B9A" if is_arc else "#B71C1C"
                svg.append(f'<line x1="{ax_px:.1f}" y1="{mg["top"]}" '
                           f'x2="{ax_px:.1f}" y2="{mg["top"]+ph}" '
                           f'stroke="{mc}" stroke-width="2.0" '
                           f'stroke-dasharray="5,3" opacity="0.9"/>')
                # diamond marker at top
                dy = mg["top"] - 6
                svg.append(f'<polygon points="{ax_px:.1f},{dy-7} {ax_px+5:.1f},{dy} '
                           f'{ax_px:.1f},{dy+7} {ax_px-5:.1f},{dy}" '
                           f'fill="{mc}"/>')
                alarm_annotations.append((ax_px, a))

    # ââ Y axis ââ
    tick_step = _nice_interval(vr, max_ticks=8)
    tv = math.ceil(vmin / tick_step) * tick_step
    while tv <= vmax + 1e-9:
        gy = ty(tv)
        if mg["top"] <= gy <= mg["top"] + ph:
            svg.append(f'<line x1="{mg["left"]-4}" y1="{gy:.1f}" '
                       f'x2="{mg["left"]}" y2="{gy:.1f}" stroke="#455A64" stroke-width="1"/>')
            svg.append(f'<text x="{mg["left"]-7}" y="{gy:.1f}" '
                       f'text-anchor="end" dominant-baseline="middle" '
                       f'font-size="10" fill="#455A64">{_fmt_v(tv)}</text>')
        tv = round(tv + tick_step, 10)

    # Y axis label
    cx = 14
    cy = mg["top"] + ph // 2
    svg.append(f'<text x="{cx}" y="{cy}" text-anchor="middle" '
               f'font-size="12" fill="#455A64" '
               f'transform="rotate(-90,{cx},{cy})">{_esc(sensor)}</text>')

    # ââ X axis ââ
    svg.append(f'<text x="{mg["left"] + pw//2}" y="{H - LABEL_H//4}" '
               f'text-anchor="middle" font-size="12" fill="#455A64">Time (s)</text>')
    tick_step_t = _nice_interval(tr, max_ticks=12)
    tt = math.ceil(t0 / tick_step_t) * tick_step_t
    while tt <= tmax + 1e-9:
        gx = tx(tt)
        if mg["left"] <= gx <= mg["left"] + pw:
            svg.append(f'<line x1="{gx:.1f}" y1="{mg["top"]+ph}" '
                       f'x2="{gx:.1f}" y2="{mg["top"]+ph+5}" stroke="#455A64" stroke-width="1"/>')
            svg.append(f'<text x="{gx:.1f}" y="{mg["top"]+ph+16}" '
                       f'text-anchor="middle" font-size="10" fill="#455A64">'
                       f'{_fmt_t(tt - t0)}</text>')
        tt = round(tt + tick_step_t, 10)

    # ââ step labels strip below X-axis ââ
    label_y0 = mg["top"] + ph + 28
    for bi, b in enumerate(clean):
        bx  = tx(b["elapsed_s"])
        ex  = tx(clean[bi + 1]["elapsed_s"]) if bi + 1 < len(clean) else tx(tmax)
        cx2 = (bx + ex) / 2
        lbl = b["step_name"][:18] if b["step_name"] else str(b["step_number"])
        svg.append(f'<text x="{cx2:.1f}" y="{label_y0}" '
                   f'text-anchor="middle" font-size="9" fill="#546E7A" '
                   f'transform="rotate(-35,{cx2:.1f},{label_y0})">'
                   f'{_esc(lbl)}</text>')

    # ââ legend (right panel) ââ
    lx = mg["left"] + pw + 18
    ly = mg["top"] + 10
    row_off = 0

    def leg_row(off, col, lw, dash, label, lcol="#37474F", bold=False):
        y = ly + off
        da = f'stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<line x1="{lx}" y1="{y+6}" x2="{lx+24}" y2="{y+6}" '
                   f'stroke="{col}" stroke-width="{lw}" {da}/>')
        fw = "bold" if bold else "normal"
        svg.append(f'<text x="{lx+28}" y="{y+10}" font-size="10" '
                   f'fill="{lcol}" font-weight="{fw}">{_esc(label)}</text>')

    # envelope legend
    svg.append(f'<rect x="{lx}" y="{ly+row_off}" width="24" height="12" '
               f'fill="#A5D6A7" opacity="0.5" stroke="#2E7D32" stroke-width="1"/>')
    svg.append(f'<text x="{lx+28}" y="{ly+row_off+10}" font-size="10" '
               f'fill="#2E7D32">Median Â± 3ÃMAD</text>')
    row_off += 18
    leg_row(row_off, "#1B5E20", "1.5", "6,3", "Median (normal)", "#1B5E20")
    row_off += 18

    # normal run lines
    for ri, rid in enumerate(sorted(normal_ids)):
        col = run_colors.get(rid, "#388E3C")
        leg_row(row_off, col, "1.2", "", f"Run {rid}  Normal", "#37474F")
        row_off += 14

    # aborted run lines
    for rid in sorted(_aborted_set_env & set(traces)):
        leg_row(row_off, "#FF6D00", "2", "", f"Run {rid}  ABORTED", "#FF6D00", bold=False)
        row_off += 14

    # anomalous run lines
    for rid in sorted(anom_ids):
        leg_row(row_off, "#D32F2F", "3.0", "", f"Run {rid}  ANOMALOUS", "#C62828", bold=True)
        row_off += 14

    # alarms / ARC events in legend
    if alarm_annotations:
        arc_anns  = [(px, a) for px, a in alarm_annotations
                     if a.get("is_arc_onset") or a.get("alarm_code") == "ARC_EVENT"]
        reg_anns  = [(px, a) for px, a in alarm_annotations
                     if not (a.get("is_arc_onset") or a.get("alarm_code") == "ARC_EVENT")]
        if arc_anns:
            row_off += 6
            leg_row(row_off, "#6A1B9A", "2.0", "5,3", "ARC onset / escalation", "#6A1B9A")
            row_off += 14
            for _, a in arc_anns:
                arc_col_s = a.get("arc_col", a["event_name"].split(".")[-1])
                level     = a.get("arc_level", "?")
                label     = f'count={level} -- {arc_col_s[:30]}'
                svg.append(f'<text x="{lx}" y="{ly+row_off+10}" font-size="9" '
                           f'fill="#6A1B9A" font-weight="bold">{_esc(label)}</text>')
                row_off += 13
        if reg_anns:
            row_off += 6
            leg_row(row_off, "#B71C1C", "1.5", "3,3", "Alarm", "#B71C1C")
            row_off += 14
            for _, a in reg_anns:
                short = a["event_name"].split(".")[-1][:28]
                svg.append(f'<text x="{lx}" y="{ly+row_off+10}" font-size="9" '
                           f'fill="#B71C1C" font-weight="bold">{_esc(short)}</text>')
                row_off += 13

    svg.append("</svg>")
    return "\n".join(svg)


# ââ B. Z-score overlay chart ââââââââââââââââââââââââââââââââââââââââââââââââââ

def svg_zscore_chart(run_id, rrows, baselines, lead, chain, tool_id,
                     alarms=None):
    """
    Z-score trace for each sensor in the lead/lag chain.
    - Shaded bands: light red above +3Ï, light blue below -3Ï
    - Lead = thick red solid line
    - Lag sensors = progressively warmer colors, thinner lines
    - Divergence onset: labeled vertical line with sensor name box
    - Alarm: bold dark-red vertical line with annotation box
    - Legend: right panel with full sensor names and divergence times
    """
    W  = ZSCORE_WIDTH
    H  = ZSCORE_HEIGHT
    LEG_W = 260
    mg = {"top": 70, "right": LEG_W + 15, "bottom": 60, "left": 85}
    pw = W - mg["left"] - mg["right"]
    ph = H - mg["top"]  - mg["bottom"]

    # build z-score traces
    plot_data = []  # (sensor, diverge_t, ts, zs, color, is_lead)
    n = len(chain)
    for ci, (sensor, div_t) in enumerate(chain):
        bl = baselines.get(sensor)
        if not bl or not bl[0]:
            continue
        ts, zs = z_trace(rrows, bl, sensor)
        if not ts:
            continue
        color = LEAD_COLOR if sensor == lead else _cascade_color(ci, n)
        plot_data.append((sensor, div_t, ts, zs, color, sensor == lead))

    if not plot_data:
        return ""

    all_t  = [t for _, _, ts, _, _, _ in plot_data for t in ts]
    all_z  = [z for _, _, _, zs, _, _ in plot_data for z in zs]
    t0, tmax = min(all_t), max(all_t)
    zpad   = max(abs(min(all_z)), abs(max(all_z)), SENSOR_DIVERGE_SIGMA + 1.5) * 1.1
    zmin, zmax = -zpad, zpad
    tr = max(tmax - t0, 1e-9)
    zr = zmax - zmin

    def tx(t): return mg["left"] + (t - t0) / tr * pw
    def ty(z): return mg["top"]  + (1 - (z - zmin) / zr) * ph

    svg = []
    svg.append(f'<svg width="{W}" height="{H}" xmlns="http://www.w3.org/2000/svg" '
               f'style="background:#FAFAFA;font-family:Arial,sans-serif;display:block;margin:10px 0">')

    # ââ title ââ
    svg.append(f'<text x="{mg["left"]+pw//2}" y="24" text-anchor="middle" '
               f'font-size="15" font-weight="bold" fill="#1A237E">'
               f'Tool {_esc(tool_id)} -- Run {_esc(run_id)} -- Z-Score Lead / Lag Cascade</text>')
    svg.append(f'<text x="{mg["left"]+pw//2}" y="42" text-anchor="middle" '
               f'font-size="11" fill="#555">'
               f'LEAD sensor: {_esc(lead)}   |   '
               f'Threshold: Â±{SENSOR_DIVERGE_SIGMA}Ï   |   '
               f'{len(plot_data)} sensors shown</text>')

    # ââ plot background ââ
    svg.append(f'<rect x="{mg["left"]}" y="{mg["top"]}" width="{pw}" height="{ph}" '
               f'fill="white" stroke="#B0BEC5" stroke-width="1.5"/>')

    # ââ grid lines ââ
    for i in range(9):
        gv = zmin + i * zr / 8
        gy = ty(gv)
        svg.append(f'<line x1="{mg["left"]}" y1="{gy:.1f}" '
                   f'x2="{mg["left"]+pw}" y2="{gy:.1f}" '
                   f'stroke="#ECEFF1" stroke-width="1"/>')

    # ââ shaded alarm bands (above +sigma and below -sigma) ââ
    y_plus  = ty(SENSOR_DIVERGE_SIGMA)
    y_minus = ty(-SENSOR_DIVERGE_SIGMA)
    y_top   = mg["top"]
    y_bot   = mg["top"] + ph
    # above +3Ï = light red
    svg.append(f'<rect x="{mg["left"]}" y="{y_top}" '
               f'width="{pw}" height="{y_plus-y_top:.1f}" '
               f'fill="#FFCDD2" opacity="0.45"/>')
    # below -3Ï = light blue
    svg.append(f'<rect x="{mg["left"]}" y="{y_minus:.1f}" '
               f'width="{pw}" height="{y_bot-y_minus:.1f}" '
               f'fill="#BBDEFB" opacity="0.45"/>')

    # ââ threshold lines ââ
    for sign in [1, -1]:
        yt = ty(sign * SENSOR_DIVERGE_SIGMA)
        color_t = "#C62828" if sign > 0 else "#1565C0"
        svg.append(f'<line x1="{mg["left"]}" y1="{yt:.1f}" '
                   f'x2="{mg["left"]+pw}" y2="{yt:.1f}" '
                   f'stroke="{color_t}" stroke-width="1.5" '
                   f'stroke-dasharray="8,5"/>')
        label_t = f'+{SENSOR_DIVERGE_SIGMA:.0f}Ï (alarm)' if sign > 0 else f'-{SENSOR_DIVERGE_SIGMA:.0f}Ï'
        svg.append(f'<text x="{mg["left"]+pw+4}" y="{yt+4:.1f}" '
                   f'font-size="10" fill="{color_t}" font-weight="bold">'
                   f'{_esc(label_t)}</text>')

    # ââ zero line ââ
    y0 = ty(0)
    svg.append(f'<line x1="{mg["left"]}" y1="{y0:.1f}" '
               f'x2="{mg["left"]+pw}" y2="{y0:.1f}" '
               f'stroke="#546E7A" stroke-width="1.2"/>')
    svg.append(f'<text x="{mg["left"]+pw+4}" y="{y0+4:.1f}" '
               f'font-size="10" fill="#546E7A">0</text>')

    # ââ Y-axis ticks ââ
    for i in range(9):
        zv = zmin + i * zr / 8
        yg = ty(zv)
        svg.append(f'<line x1="{mg["left"]-6}" y1="{yg:.1f}" '
                   f'x2="{mg["left"]}" y2="{yg:.1f}" '
                   f'stroke="#90A4AE" stroke-width="1.5"/>')
        svg.append(f'<text x="{mg["left"]-10}" y="{yg+4:.1f}" '
                   f'text-anchor="end" font-size="11" fill="#37474F">{zv:.1f}</text>')
    midy = mg["top"] + ph / 2
    svg.append(f'<text x="16" y="{midy:.1f}" text-anchor="middle" '
               f'font-size="12" fill="#455A64" '
               f'transform="rotate(-90,16,{midy:.1f})">'
               f'Z-Score (sigma from baseline)</text>')

    # ââ X-axis ticks ââ
    ti = _nice_interval(tr)
    tt = 0.0
    while tt <= tmax + 0.5:
        xp = tx(tt)
        svg.append(f'<line x1="{xp:.1f}" y1="{mg["top"]+ph}" '
                   f'x2="{xp:.1f}" y2="{mg["top"]+ph+6}" '
                   f'stroke="#90A4AE" stroke-width="1.5"/>')
        svg.append(f'<text x="{xp:.1f}" y="{mg["top"]+ph+20}" '
                   f'text-anchor="middle" font-size="11" fill="#37474F">'
                   f'{_fmt_t(tt)}</text>')
        tt += ti
    svg.append(f'<text x="{mg["left"]+pw//2}" y="{H-8}" '
               f'text-anchor="middle" font-size="12" fill="#455A64" font-weight="bold">'
               f'Recipe Elapsed Time (s)</text>')

    # ââ traces: lags first (behind), lead on top ââ
    for sensor, div_t, ts, zs, color, is_lead in sorted(plot_data, key=lambda x: x[5]):
        if is_lead:
            continue  # draw last
        pts = " ".join(f"{tx(t):.1f},{ty(z):.1f}" for t, z in zip(ts, zs))
        svg.append(f'<polyline points="{pts}" fill="none" '
                   f'stroke="{color}" stroke-width="2" opacity="0.80"/>')

    # lead trace: white halo + red line
    for sensor, div_t, ts, zs, color, is_lead in plot_data:
        if not is_lead:
            continue
        pts = " ".join(f"{tx(t):.1f},{ty(z):.1f}" for t, z in zip(ts, zs))
        svg.append(f'<polyline points="{pts}" fill="none" '
                   f'stroke="white" stroke-width="7" opacity="0.6"/>')
        svg.append(f'<polyline points="{pts}" fill="none" '
                   f'stroke="{LEAD_COLOR}" stroke-width="3.5" opacity="1.0"/>')

    # ââ divergence onset markers with labeled boxes ââ
    for ci, (sensor, div_t, ts, zs, color, is_lead) in enumerate(plot_data):
        xd = tx(div_t)
        svg.append(f'<line x1="{xd:.1f}" y1="{mg["top"]}" '
                   f'x2="{xd:.1f}" y2="{mg["top"]+ph}" '
                   f'stroke="{color}" stroke-width="1.5" '
                   f'stroke-dasharray="4,4" opacity="0.75"/>')
        # labeled box at top, staggered vertically to avoid overlap
        box_y = mg["top"] + 4 + ci * 20
        tag   = "LEAD" if is_lead else f"lag-{ci}"
        lbl   = f'{tag}: {div_t:.0f}s'
        bw    = len(lbl) * 6 + 8
        svg.append(f'<rect x="{xd+2:.1f}" y="{box_y}" '
                   f'width="{bw}" height="14" rx="3" '
                   f'fill="{color}" opacity="0.9"/>')
        svg.append(f'<text x="{xd+6:.1f}" y="{box_y+10}" '
                   f'font-size="9" fill="white" font-weight="bold">'
                   f'{_esc(lbl)}</text>')

    # ââ alarm markers ââ
    if alarms:
        for a in alarms:
            at = a["elapsed_s"]
            if at < t0 or at > tmax:
                continue
            ax = tx(at)
            svg.append(f'<line x1="{ax:.1f}" y1="{mg["top"]}" '
                       f'x2="{ax:.1f}" y2="{mg["top"]+ph}" '
                       f'stroke="#B71C1C" stroke-width="3" '
                       f'stroke-dasharray="10,4" opacity="0.95"/>')
            # annotation box at bottom
            short = a["event_name"].split(".")[-1][:28]
            desc  = a["description"][:45].strip()
            bw2   = max(len(short), len(desc)) * 6 + 10
            bx2   = min(ax + 4, mg["left"] + pw - bw2 - 4)
            by2   = mg["top"] + ph - 52
            svg.append(f'<rect x="{bx2:.1f}" y="{by2}" '
                       f'width="{bw2}" height="46" rx="4" '
                       f'fill="#FFEBEE" stroke="#B71C1C" stroke-width="1.5"/>')
            svg.append(f'<text x="{bx2+6:.1f}" y="{by2+14}" '
                       f'font-size="10" fill="#B71C1C" font-weight="bold">'
                       f'FAULT ALARM</text>')
            svg.append(f'<text x="{bx2+6:.1f}" y="{by2+28}" '
                       f'font-size="9" fill="#B71C1C">{_esc(short)}</text>')
            svg.append(f'<text x="{bx2+6:.1f}" y="{by2+40}" '
                       f'font-size="9" fill="#555">{_esc(desc)}</text>')

    # ââ right legend panel ââ
    lx = mg["left"] + pw + 15
    ly = mg["top"]
    row_y = ly

    def z_leg(color, lw, label, sublabel="", bold=False):
        nonlocal row_y
        svg.append(f'<line x1="{lx}" y1="{row_y+8}" x2="{lx+30}" y2="{row_y+8}" '
                   f'stroke="{color}" stroke-width="{lw}"/>')
        fw2 = "bold" if bold else "normal"
        svg.append(f'<text x="{lx+36}" y="{row_y+12}" font-size="11" '
                   f'fill="{color}" font-weight="{fw2}">{_esc(label)}</text>')
        if sublabel:
            svg.append(f'<text x="{lx+36}" y="{row_y+24}" font-size="9" '
                       f'fill="#888">{_esc(sublabel)}</text>')
        row_y += 28 if sublabel else 20

    row_y += 4
    # shaded band legend
    svg.append(f'<rect x="{lx}" y="{row_y}" width="30" height="14" '
               f'fill="#FFCDD2" opacity="0.7" stroke="#C62828" stroke-width="1"/>')
    svg.append(f'<text x="{lx+36}" y="{row_y+11}" font-size="11" fill="#C62828">'
               f'Above +{SENSOR_DIVERGE_SIGMA:.0f}Ï (alarm zone)</text>')
    row_y += 22
    svg.append(f'<rect x="{lx}" y="{row_y}" width="30" height="14" '
               f'fill="#BBDEFB" opacity="0.7" stroke="#1565C0" stroke-width="1"/>')
    svg.append(f'<text x="{lx+36}" y="{row_y+11}" font-size="11" fill="#1565C0">'
               f'Below -{SENSOR_DIVERGE_SIGMA:.0f}Ï</text>')
    row_y += 28

    svg.append(f'<line x1="{lx}" y1="{row_y}" x2="{lx+30}" y2="{row_y}" '
               f'stroke="#546E7A" stroke-width="1.2"/>')
    svg.append(f'<text x="{lx+36}" y="{row_y+4}" font-size="11" fill="#546E7A">Zero (baseline)</text>')
    row_y += 22

    # sensor entries
    for ci, (sensor, div_t, ts, zs, color, is_lead) in enumerate(plot_data):
        lw   = "3.5" if is_lead else "2"
        tag  = "LEAD" if is_lead else f"lag-{ci}"
        z_leg(color, lw, f"[{tag}]  {sensor}",
              f"diverges at {div_t:.1f}s", bold=is_lead)

    # alarm legend entry
    if alarms:
        row_y += 6
        svg.append(f'<line x1="{lx}" y1="{row_y+8}" x2="{lx+30}" y2="{row_y+8}" '
                   f'stroke="#B71C1C" stroke-width="3" stroke-dasharray="8,3"/>')
        svg.append(f'<text x="{lx+36}" y="{row_y+12}" font-size="11" '
                   f'fill="#B71C1C" font-weight="bold">Fault alarm event</text>')

    svg.append("</svg>")
    return "\n".join(svg)


# ââ C. SP deviation chart: actual vs setpoint for one sensor pair ââââââââââââ

def svg_sp_deviation(actual_col, sp_col, run_dict, sp_stats,
                     run_colors, anomalous_runs):
    """
    Overlay chart showing actual sensor value and its setpoint for every run.
    - Setpoint (SP): drawn as a DASHED line in the same color as the run
    - Actual value:  drawn as a SOLID line in the same color
    - Off-SP region: shaded in light red for anomalous runs
    - Anomalous run actual: thick red solid line (same as raw trace convention)
    - Normal envelope (actual): light green band
    - Secondary Y-axis band shows % deviation from SP
    """
    # collect traces
    act_traces = {}    # rid -> (ts, vs)
    sp_traces  = {}    # rid -> (ts, vs)
    for rid, rrows in run_dict.items():
        sp_d = sp_stats.get(rid, {}).get((actual_col, sp_col))
        if not sp_d or sp_d["n_active"] < 5:
            continue
        act_traces[rid] = (sp_d["act_ts"], sp_d["act_vs"])
        sp_traces[rid]  = (sp_d["sp_ts"],  sp_d["sp_vs"])

    if not act_traces:
        return ""

    # value bounds (both actual and SP)
    all_ts = [t for ts, _ in act_traces.values() for t in ts]
    all_vs = ([v for _, vs in act_traces.values() for v in vs] +
              [v for _, vs in sp_traces.values()  for v in vs])
    t0, tmax = min(all_ts), max(all_ts)
    vmin_d, vmax_d = min(all_vs), max(all_vs)
    pad  = max((vmax_d - vmin_d) * 0.12, abs(vmax_d) * 0.05, 0.1)
    vmin = vmin_d - pad;  vmax = vmax_d + pad
    tr   = max(tmax - t0, 1e-9)
    vr   = max(vmax - vmin, 1e-9)

    LEG_W = 230
    mg = {"top": 65, "right": LEG_W + 20, "bottom": 55, "left": 85}
    W  = max(900, min(MAX_CHART_WIDTH, 1400))
    H  = 360
    pw = W - mg["left"] - mg["right"]
    ph = H - mg["top"]  - mg["bottom"]

    def tx(t): return mg["left"] + (t - t0) / tr * pw
    def ty(v): return mg["top"]  + (1 - (v - vmin) / vr) * ph

    svg = []
    svg.append(f'<svg width="{W}" height="{H}" xmlns="http://www.w3.org/2000/svg" '
               f'style="background:#FAFAFA;font-family:Arial,sans-serif;display:block;margin:6px 0">')

    # title
    svg.append(f'<text x="{mg["left"]+pw//2}" y="22" text-anchor="middle" '
               f'font-size="14" font-weight="bold" fill="#1A237E">'
               f'{_esc(actual_col)}  vs  Setpoint ({_esc(sp_col)})</text>')
    anom_label = ", ".join(f"Run {r}" for r in sorted(anomalous_runs) if r in act_traces)
    svg.append(f'<text x="{mg["left"]+pw//2}" y="40" text-anchor="middle" '
               f'font-size="11" fill="#555">'
               f'Solid = actual  |  Dashed = setpoint  |  '
               f'Red shading = anomalous run off-SP  |  Anomalous: {_esc(anom_label or "none")}'
               f'</text>')

    # plot background
    svg.append(f'<rect x="{mg["left"]}" y="{mg["top"]}" width="{pw}" height="{ph}" '
               f'fill="white" stroke="#B0BEC5" stroke-width="1.5"/>')

    # grid
    for i in range(7):
        gv = vmin + i * vr / 6
        gy = ty(gv)
        svg.append(f'<line x1="{mg["left"]}" y1="{gy:.1f}" x2="{mg["left"]+pw}" y2="{gy:.1f}" '
                   f'stroke="#ECEFF1" stroke-width="1"/>')

    # normal envelope (actual values of normal runs)
    normal_ids = [r for r in act_traces if r not in anomalous_runs]
    if len(normal_ids) >= 2:
        from collections import defaultdict as _dd
        BIN = 1.0
        env = _dd(list)
        for rid in normal_ids:
            for t, v in zip(*act_traces[rid]):
                env[round(t / BIN) * BIN].append(v)
        env_t  = sorted(env)
        env_lo = [min(env[t]) for t in env_t]
        env_hi = [max(env[t]) for t in env_t]
        upper  = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(env_t, env_hi))
        lower  = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(reversed(env_t), reversed(env_lo)))
        svg.append(f'<polygon points="{upper} {lower}" fill="#A5D6A7" opacity="0.30"/>')

    # off-SP shading for anomalous runs
    for rid in sorted(anomalous_runs):
        sp_d = sp_stats.get(rid, {}).get((actual_col, sp_col))
        if not sp_d or not sp_d["off_ts"]:
            continue
        # shade regions where actual is off-SP: fill between actual and SP line
        # Build combined list, sorted by time
        pts_above, pts_below = [], []
        for t, v_act, v_sp in zip(sp_d["act_ts"], sp_d["act_vs"], sp_d["sp_vs"]):
            if abs(v_act - v_sp) / max(abs(v_sp), 1e-9) > SP_OFF_FRAC_THRESHOLD:
                if v_act > v_sp:
                    pts_above.append((t, v_act, v_sp))
                else:
                    pts_below.append((t, v_act, v_sp))
        for pts, color in [(pts_above, "#FFCDD2"), (pts_below, "#BBDEFB")]:
            if len(pts) < 2:
                continue
            top_pts = " ".join(f"{tx(t):.1f},{ty(va):.1f}" for t, va, _ in pts)
            bot_pts = " ".join(f"{tx(t):.1f},{ty(vs):.1f}" for t, _, vs in reversed(pts))
            svg.append(f'<polygon points="{top_pts} {bot_pts}" '
                       f'fill="{color}" opacity="0.55"/>')

    # SP lines (dashed, same color as run, thinner)
    for rid, (ts, vs) in sorted(sp_traces.items()):
        col = run_colors.get(rid, "#777")
        pts = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(ts, vs))
        lw  = "2" if rid not in anomalous_runs else "2"
        svg.append(f'<polyline points="{pts}" fill="none" stroke="{col}" '
                   f'stroke-width="{lw}" stroke-dasharray="10,5" opacity="0.70"/>')

    # actual lines -- normal runs thin, anomalous thick red with halo
    for rid in sorted(act_traces):
        if rid in anomalous_runs:
            continue
        ts, vs = act_traces[rid]
        col = run_colors.get(rid, "#777")
        pts = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(ts, vs))
        svg.append(f'<polyline points="{pts}" fill="none" stroke="{col}" '
                   f'stroke-width="1.8" opacity="0.65"/>')

    for rid in sorted(anomalous_runs):
        if rid not in act_traces:
            continue
        ts, vs = act_traces[rid]
        pts = " ".join(f"{tx(t):.1f},{ty(v):.1f}" for t, v in zip(ts, vs))
        svg.append(f'<polyline points="{pts}" fill="none" stroke="white" '
                   f'stroke-width="6" opacity="0.65"/>')
        svg.append(f'<polyline points="{pts}" fill="none" stroke="#D32F2F" '
                   f'stroke-width="3" opacity="1.0"/>')

    # Y axis
    for i in range(7):
        v  = vmin + i * vr / 6
        yp = ty(v)
        svg.append(f'<line x1="{mg["left"]-6}" y1="{yp:.1f}" x2="{mg["left"]}" y2="{yp:.1f}" '
                   f'stroke="#90A4AE" stroke-width="1.5"/>')
        svg.append(f'<text x="{mg["left"]-10}" y="{yp+4:.1f}" text-anchor="end" '
                   f'font-size="11" fill="#37474F">{_fmt_v(v)}</text>')
    midy = mg["top"] + ph / 2
    svg.append(f'<text x="16" y="{midy:.1f}" text-anchor="middle" font-size="11" fill="#455A64" '
               f'transform="rotate(-90,16,{midy:.1f})">{_esc(actual_col)}</text>')

    # X axis
    ti = _nice_interval(tr)
    tt = 0.0
    while tt <= tmax + 0.5:
        xp = tx(tt)
        svg.append(f'<line x1="{xp:.1f}" y1="{mg["top"]+ph}" x2="{xp:.1f}" y2="{mg["top"]+ph+6}" '
                   f'stroke="#90A4AE" stroke-width="1.5"/>')
        svg.append(f'<text x="{xp:.1f}" y="{mg["top"]+ph+20}" text-anchor="middle" '
                   f'font-size="11" fill="#37474F">{_fmt_t(tt)}</text>')
        tt += ti
    svg.append(f'<text x="{mg["left"]+pw//2}" y="{H-5}" text-anchor="middle" '
               f'font-size="11" fill="#455A64" font-weight="bold">Recipe Elapsed Time (s)</text>')

    # legend
    lx = mg["left"] + pw + 18;  row_y = mg["top"]
    svg.append(f'<rect x="{lx}" y="{row_y}" width="12" height="12" '
               f'fill="#A5D6A7" opacity="0.5" stroke="#388E3C" stroke-width="1"/>')
    svg.append(f'<text x="{lx+18}" y="{row_y+11}" font-size="10" fill="#2E7D32">'
               f'Normal envelope</text>')
    row_y += 20
    svg.append(f'<rect x="{lx}" y="{row_y}" width="12" height="12" '
               f'fill="#FFCDD2" opacity="0.7"/>')
    svg.append(f'<text x="{lx+18}" y="{row_y+11}" font-size="10" fill="#C62828">'
               f'Actual above SP (&gt;5%)</text>')
    row_y += 20
    svg.append(f'<rect x="{lx}" y="{row_y}" width="12" height="12" '
               f'fill="#BBDEFB" opacity="0.7"/>')
    svg.append(f'<text x="{lx+18}" y="{row_y+11}" font-size="10" fill="#1565C0">'
               f'Actual below SP (&gt;5%)</text>')
    row_y += 24

    for rid in sorted(act_traces):
        col   = run_colors.get(rid, "#777")
        is_a  = rid in anomalous_runs
        lw_act = "3" if is_a else "1.8"
        lw_sp  = "2"
        clr    = "#D32F2F" if is_a else col
        # actual line swatch
        svg.append(f'<line x1="{lx}" y1="{row_y+7}" x2="{lx+22}" y2="{row_y+7}" '
                   f'stroke="{clr}" stroke-width="{lw_act}"/>')
        # SP dashed swatch
        svg.append(f'<line x1="{lx}" y1="{row_y+15}" x2="{lx+22}" y2="{row_y+15}" '
                   f'stroke="{col}" stroke-width="{lw_sp}" stroke-dasharray="6,3"/>')
        label = f'Run {rid}{"  ANOM" if is_a else ""}'
        svg.append(f'<text x="{lx+28}" y="{row_y+10}" font-size="10" '
                   f'fill="{"#C62828" if is_a else col}" '
                   f'font-weight="{"bold" if is_a else "normal"}">{_esc(label)}</text>')
        svg.append(f'<text x="{lx+28}" y="{row_y+22}" font-size="9" fill="#888">'
                   f'âââ actual   - - - SP</text>')
        row_y += 30

    svg.append("</svg>")
    return "\n".join(svg)


# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
# 7. HTML REPORT
# ââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ

def svg_trend_chart(sensor, step_num, run_means_for_sensor,
                    run_order, run_tool, tools,
                    aborted_runs=None, row_counts=None,
                    highflier_info=None):
    """
    Line+scatter chart: x = run index, y = per-run mean for this sensor.
    Aborted runs shown as orange triangles at the bottom of the chart area.
    Highfliers shown as coloured diamonds (red = real, amber = artifact).
    """
    W = TREND_CHART_W; H = TREND_CHART_H
    mt = TREND_MARGIN['top']; mr = TREND_MARGIN['right']
    mb = TREND_MARGIN['bottom']; ml = TREND_MARGIN['left']
    pw = W - ml - mr; ph = H - mt - mb

    aborted_runs  = aborted_runs  or set()
    row_counts    = row_counts    or {}
    highflier_info = highflier_info or {}

    tool_color = {t: TOOL_COLORS[i % len(TOOL_COLORS)]
                  for i, t in enumerate(sorted(tools))}

    # Only non-aborted runs that have data (run_order contains composite keys)
    pts = [(i, ckey, run_tool.get(ckey, 'ALL'), run_means_for_sensor[ckey])
           for i, ckey in enumerate(run_order)
           if ckey in run_means_for_sensor and ckey not in aborted_runs]
    if not pts:
        return ''

    vals     = [v for _,_,_,v in pts]
    ymin     = min(vals); ymax = max(vals)
    ypad     = (ymax - ymin) * 0.15 if ymax > ymin else max(abs(ymax)*0.1, 1.0)
    ylo      = ymin - ypad; yhi = ymax + ypad

    n_runs   = len(run_order)
    def xp(i): return ml + pw * i / max(n_runs - 1, 1)
    def yp(v): return mt + ph * (1 - (v - ylo) / max(yhi - ylo, 1e-12))

    # Per-tool statistics for separate bands
    tool_stats = {}   # tool -> (mean, std)
    by_tool_vals = defaultdict(list)
    for _, ckey, tool, v in pts:
        by_tool_vals[tool].append(v)
    for tool, tvals in by_tool_vals.items():
        if len(tvals) >= 2:
            tm = sum(tvals) / len(tvals)
            ts = math.sqrt(sum((v - tm) ** 2 for v in tvals) / len(tvals))
            tool_stats[tool] = (tm, ts)

    mean_all = sum(vals) / len(vals)
    std_all  = math.sqrt(sum((v-mean_all)**2 for v in vals) / len(vals))

    svg = [f'<svg width="{W}" height="{H}" '
           f'style="font-family:Arial,sans-serif;display:block;margin-bottom:4px">']

    # plot background
    svg.append(f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" '
               f'fill="#FAFAFA" stroke="#DDD" stroke-width="1"/>')

    # y-axis grid + labels
    def _nice_ticks(lo, hi, n=5):
        span = hi - lo
        if span <= 0: return [lo]
        raw = span / n
        mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
        for f in (1, 2, 2.5, 5, 10):
            step = f * mag
            start = math.ceil(lo / step) * step
            ticks = []; t = start
            while t <= hi + step * 0.01:
                ticks.append(round(t, 10)); t += step
            if 3 <= len(ticks) <= 8:
                return ticks
        return [lo, hi]

    for t in _nice_ticks(ylo, yhi):
        yy = yp(t)
        if mt <= yy <= mt + ph:
            svg.append(f'<line x1="{ml}" y1="{yy:.1f}" x2="{ml+pw}" y2="{yy:.1f}" '
                       f'stroke="#E0E0E0" stroke-width="1"/>')
            lbl = f'{t:.1f}' if abs(t) < 1000 else f'{t:.0f}'
            svg.append(f'<text x="{ml-5}" y="{yy+4:.1f}" text-anchor="end" '
                       f'font-size="10" fill="#666">{lbl}</text>')

    # Per-tool mean +/- std bands (each tool gets its own color band)
    band_fill_opacity = 0.12 if len(tool_stats) > 1 else 0.30
    for tool in sorted(tool_stats):
        tm, ts = tool_stats[tool]
        tcol = tool_color.get(tool, '#1565C0')
        b1 = yp(tm + ts); b2 = yp(tm - ts)
        svg.append(f'<rect x="{ml}" y="{min(b1,b2):.1f}" width="{pw}" '
                   f'height="{abs(b2-b1):.1f}" fill="{tcol}" opacity="{band_fill_opacity}"/>')
        ym = yp(tm)
        svg.append(f'<line x1="{ml}" y1="{ym:.1f}" x2="{ml+pw}" y2="{ym:.1f}" '
                   f'stroke="{tcol}" stroke-width="1.2" stroke-dasharray="6,3" opacity="0.7"/>')
    # per-tool connecting lines (non-aborted only)
    by_tool = defaultdict(list)
    for i, ckey, tool, v in pts:
        by_tool[tool].append((i, v))
    for tool, tpts in by_tool.items():
        col = tool_color.get(tool, '#555')
        tpts_s = sorted(tpts)
        if len(tpts_s) >= 2:
            coords = ' '.join(f'{xp(i):.1f},{yp(v):.1f}' for i, v in tpts_s)
            svg.append(f'<polyline points="{coords}" fill="none" '
                       f'stroke="{col}" stroke-width="1.4" opacity="0.45"/>')

    # normal dots + highflier diamonds
    for i, ckey, tool, v in pts:
        rid_lbl = ckey.split("::", 1)[1] if "::" in ckey else ckey
        col = tool_color.get(tool, '#555')
        hf  = highflier_info.get(ckey)
        if hf:
            # diamond shape
            cx = xp(i); cy = yp(v); r = 5.5
            hf_col = HIGHFLIER_REAL if hf['class'] == 'real' else HIGHFLIER_ART
            pts_d  = (f'{cx},{cy-r} {cx+r},{cy} {cx},{cy+r} {cx-r},{cy}')
            tip = (f'Run {rid_lbl} HIGHFLIER ({hf["direction"]}) '
                   f'val={v:.4g} | {hf["class"].upper()}: {hf["reason"]}')
            svg.append(f'<polygon points="{pts_d}" fill="{hf_col}" '
                       f'stroke="white" stroke-width="1">'
                       f'<title>{_esc(tip)}</title></polygon>')
        else:
            svg.append(f'<circle cx="{xp(i):.1f}" cy="{yp(v):.1f}" r="3.5" '
                       f'fill="{col}" stroke="white" stroke-width="0.8" opacity="0.85">'
                       f'<title>Run {_esc(rid_lbl)} ({_esc(tool)}): {v:.4g}</title>'
                       f'</circle>')

    # aborted run markers -- orange triangles along bottom edge
    abort_y = mt + ph - 8
    for i, ckey in enumerate(run_order):
        if ckey not in aborted_runs:
            continue
        xx      = xp(i)
        rid_lbl = ckey.split("::", 1)[1] if "::" in ckey else ckey
        cnt     = row_counts.get(ckey, '?')
        svg.append(f'<polygon points="{xx},{abort_y} {xx-5},{abort_y+12} '
                   f'{xx+5},{abort_y+12}" fill="{ABORT_COLOR}" opacity="0.9">'
                   f'<title>Run {_esc(rid_lbl)} ABORTED ({cnt} rows)</title>'
                   f'</polygon>')
        svg.append(f'<text x="{xx}" y="{abort_y+22}" text-anchor="middle" '
                   f'font-size="8" fill="{ABORT_COLOR}">ABORT</text>')

    # axes
    svg.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" '
               f'stroke="#999" stroke-width="1.5"/>')
    svg.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" '
               f'stroke="#999" stroke-width="1.5"/>')

    # x-axis ticks (show only the rid part, strip tool:: prefix)
    tick_every = max(1, n_runs // 10)
    for i, ckey in enumerate(run_order):
        if i % tick_every == 0 or i == n_runs - 1:
            xx   = xp(i)
            lbl  = ckey.split("::", 1)[1] if "::" in ckey else ckey
            svg.append(f'<line x1="{xx:.1f}" y1="{mt+ph}" x2="{xx:.1f}" '
                       f'y2="{mt+ph+5}" stroke="#999" stroke-width="1"/>')
            svg.append(f'<text x="{xx:.1f}" y="{mt+ph+16}" text-anchor="middle" '
                       f'font-size="9" fill="#666">{_esc(lbl)}</text>')

    # axis labels
    svg.append(f'<text x="{ml+pw//2}" y="{H-6}" text-anchor="middle" '
               f'font-size="10" fill="#555">Run (chronological order)</text>')
    svg.append(f'<text transform="rotate(-90,{ml-48},{mt+ph//2})" '
               f'x="{ml-48}" y="{mt+ph//2+4}" text-anchor="middle" '
               f'font-size="10" fill="#555">{_esc(_short(sensor))}</text>')

    # chart title -- show per-tool mean/std when multiple tools present
    svg.append(f'<text x="{ml}" y="{mt-10}" font-size="11" font-weight="bold" '
               f'fill="#333">{_esc(sensor)}</text>')
    if len(tool_stats) > 1:
        stats_parts = []
        for _t in sorted(tool_stats):
            _tm, _ts = tool_stats[_t]
            _short_t = "_".join(_t.split("_")[-2:]) if "_" in _t else _t
            _cv = f'{_ts/abs(_tm):.3f}' if abs(_tm) > 1e-9 else 'n/a'
            stats_parts.append(f'{_short_t}: mean={_tm:.4g} std={_ts:.4g} CV={_cv}')
        stats_str = "   |   ".join(stats_parts)
    else:
        cv_str = f'{std_all/abs(mean_all):.3f}' if abs(mean_all) > 1e-9 else 'n/a'
        stats_str = f'mean={mean_all:.4g}  std={std_all:.4g}  CV={cv_str}'
    svg.append(f'<text x="{ml+pw}" y="{mt-10}" text-anchor="end" font-size="10" '
               f'fill="#888">{_esc(stats_str)}</text>')

    # legend row: tools + highflier symbols + abort symbol
    lx = ml; ly = mt + ph + 34
    for tool in sorted(tools):
        col = tool_color.get(tool, '#555')
        _parts = tool.split("_")
        short_tool = "_".join(_parts[-3:]) if len(_parts) >= 3 else tool
        # colored dot for the tool's runs
        svg.append(f'<rect x="{lx}" y="{ly}" width="11" height="8" fill="{col}"/>')
        # band swatch (semi-transparent fill matching the band)
        band_op = 0.25 if len(tool_stats) > 1 else 0.50
        svg.append(f'<rect x="{lx+13}" y="{ly}" width="22" height="8" '
                   f'fill="{col}" opacity="{band_op}" stroke="{col}" stroke-width="0.5"/>')
        svg.append(f'<text x="{lx+38}" y="{ly+8}" font-size="9" fill="#444">'
                   f'{_esc(short_tool)} (mean+/-std band)</text>')
        lx += min(200, pw // max(len(tools) + 2, 1))
    # highflier real
    svg.append(f'<polygon points="{lx},{ly+4} {lx+5},{ly} {lx+10},{ly+4} '
               f'{lx+5},{ly+8}" fill="{HIGHFLIER_REAL}"/>')
    svg.append(f'<text x="{lx+13}" y="{ly+8}" font-size="9" fill="#444">'
               f'Highflier-Real</text>')
    lx += 110
    # highflier artifact
    svg.append(f'<polygon points="{lx},{ly+4} {lx+5},{ly} {lx+10},{ly+4} '
               f'{lx+5},{ly+8}" fill="{HIGHFLIER_ART}"/>')
    svg.append(f'<text x="{lx+13}" y="{ly+8}" font-size="9" fill="#444">'
               f'Highflier-Artifact</text>')
    lx += 120
    # abort marker
    svg.append(f'<polygon points="{lx+5},{ly} {lx},{ly+8} {lx+10},{ly+8}" '
               f'fill="{ABORT_COLOR}"/>')
    svg.append(f'<text x="{lx+13}" y="{ly+8}" font-size="9" fill="#444">'
               f'Aborted run</text>')

    svg.append('</svg>')
    return '\n'.join(svg)

def svg_corr_heatmap(sensors, corr, step_name, min_abs_corr=0.5):
    """
    Render a Pearson correlation heatmap as an SVG.
    Color scale: -1 -> blue, 0 -> white, +1 -> red.

    Only sensors with at least one |r| >= min_abs_corr off-diagonal are shown.
    Cells below the threshold are rendered as blank white.
    Full sensor names are printed without truncation.
    """
    # Filter to sensors with at least one strong off-diagonal correlation
    keep_idx = [i for i in range(len(sensors))
                if any(abs(corr[i][j]) >= min_abs_corr
                       for j in range(len(sensors)) if j != i)]
    if not keep_idx:
        return ''
    sensors = [sensors[i] for i in keep_idx]
    corr    = [[corr[i][j] for j in keep_idx] for i in keep_idx]

    ns   = len(sensors)
    if ns == 0:
        return ''
    cell  = max(36, min(64, 600 // ns))
    # Compute label padding from the longest full sensor name
    max_label_len = max(len(s) for s in sensors)
    pad_l = max(200, min(600, max_label_len * 7))
    pad_t = 20    # top padding
    pad_b = max(200, min(600, max_label_len * 7))  # bottom for rotated x labels
    pad_r = 20
    W     = pad_l + ns * cell + pad_r
    H     = pad_t + ns * cell + pad_b


    def _color(r):
        r = max(-1.0, min(1.0, r))
        if r >= 0:
            w = 1.0 - r
            return f'rgb({int(255*w)},{int(255*w)},255)' if False else \
                   f'rgb({int(255-(255-180)*r)},{int(255-(255-50)*r)},{int(255-(255-50)*r)})'
        else:
            w = 1.0 + r   # 0=full blue, 1=white
            return f'rgb({int(255*w)},{int(255-(255-50)*(-r))},{int(255-(255-50)*(-r))})'

    def _color2(r):
        """Blue(-1) â White(0) â Red(+1)"""
        r = max(-1.0, min(1.0, r))
        if r >= 0:
            g = int(255 * (1 - r)); b = int(255 * (1 - r))
            return f'rgb(255,{g},{b})'
        else:
            g = int(255 * (1 + r)); r2 = int(255 * (1 + r))
            return f'rgb({r2},{g},255)'

    svg = [f'<svg width="{W}" height="{H}" '
           f'style="font-family:Arial,sans-serif;display:block">']

    for i in range(ns):
        for j in range(ns):
            x  = pad_l + j * cell
            y  = pad_t + i * cell
            r  = corr[i][j]
            # Show color only when |r| >= threshold (or diagonal); blank otherwise
            if i == j or abs(r) >= min_abs_corr:
                fc = _color2(r)
            else:
                fc = '#F5F5F5'   # near-white for below-threshold cells
            svg.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" '
                       f'fill="{fc}" stroke="white" stroke-width="0.5"/>')
            # Only print value when |r| >= threshold or diagonal
            if i == j or abs(r) >= min_abs_corr:
                txt_col = '#333' if abs(r) < 0.6 else 'white'
                lbl     = f'{r:.2f}' if i != j else '1.00'
                svg.append(f'<text x="{x+cell//2}" y="{y+cell//2+4}" '
                           f'text-anchor="middle" font-size="{max(8, cell//4+2)}" '
                           f'fill="{txt_col}">{lbl}</text>')

    # y-axis labels -- full sensor name, no truncation
    for i, s in enumerate(sensors):
        y = pad_t + i * cell + cell // 2 + 4
        svg.append(f'<text x="{pad_l-6}" y="{y}" text-anchor="end" '
                   f'font-size="11" fill="#333">{_esc(s)}</text>')

    # x-axis labels (rotated -45 deg) -- full sensor name, no truncation
    for j, s in enumerate(sensors):
        x = pad_l + j * cell + cell // 2
        y = pad_t + ns * cell + 8
        svg.append(f'<text transform="rotate(-45,{x},{y})" x="{x}" y="{y}" '
                   f'text-anchor="end" font-size="11" fill="#333">{_esc(s)}</text>')

    # color scale legend
    leg_x = pad_l; leg_y = H - 18; leg_w = ns * cell; leg_h = 10
    steps = 40
    for k in range(steps):
        r_val = -1.0 + 2.0 * k / steps
        fc    = _color2(r_val)
        bx    = leg_x + int(leg_w * k / steps)
        bw    = max(1, int(leg_w / steps) + 1)
        svg.append(f'<rect x="{bx}" y="{leg_y}" width="{bw}" height="{leg_h}" '
                   f'fill="{fc}"/>')
    svg.append(f'<text x="{leg_x}" y="{leg_y-3}" font-size="9" fill="#666">-1</text>')
    svg.append(f'<text x="{leg_x+leg_w//2}" y="{leg_y-3}" text-anchor="middle" '
               f'font-size="9" fill="#666">0</text>')
    svg.append(f'<text x="{leg_x+leg_w}" y="{leg_y-3}" text-anchor="end" '
               f'font-size="9" fill="#666">+1</text>')

    svg.append('</svg>')
    return '\n'.join(svg)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7  --  HTML REPORT
# ══════════════════════════════════════════════════════════════════════════════

HTML_STYLE = """<style>
body{font-family:Arial,sans-serif;margin:24px;background:#f4f4f4;color:#333}
h1{color:#1A237E;border-bottom:3px solid #1565C0;padding-bottom:8px;margin-bottom:6px}
h2{color:#1565C0;margin-top:32px;border-left:4px solid #1565C0;padding-left:10px}
h3{color:#555;margin-top:18px;margin-bottom:6px}
.card{background:white;border-radius:6px;padding:18px;margin:14px 0;
      box-shadow:0 2px 6px rgba(0,0,0,0.09);overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{border:1px solid #ddd;padding:7px 10px;text-align:left}
th{background:#E3F2FD;font-weight:bold}
.anom td{background:#FFEBEE}
.norm td{background:#E8F5E9}
.abort td{background:#FFF3E0}
.top1 td{background:#FFF8E1}
.real td{background:#FFEBEE}
.artifact td{background:#FFFDE7}
.chip-a{display:inline-block;background:#C62828;color:white;border-radius:3px;padding:1px 7px;font-size:11px}
.chip-n{display:inline-block;background:#2E7D32;color:white;border-radius:3px;padding:1px 7px;font-size:11px}
.chip-ab{display:inline-block;background:#E65100;color:white;border-radius:3px;padding:1px 7px;font-size:11px}
.badge-top{display:inline-block;background:#1565C0;color:white;border-radius:3px;padding:1px 7px;font-size:11px;margin-left:6px}
.feat-grid{display:grid;grid-template-columns:1fr 1fr;gap:12px}
.feat-card{background:#F8F9FA;border:1px solid #DEE2E6;border-radius:5px;padding:14px}
.feat-card h4{color:#1A237E;margin:0 0 6px 0;font-size:13px}
.feat-card p{margin:0;font-size:12px;color:#555;line-height:1.5}
.sensor-grid{display:grid;grid-template-columns:1fr 1fr;gap:14px}
svg{display:block}
details summary{cursor:pointer;font-weight:bold;color:#1565C0;padding:6px 0}
details[open] summary{color:#C62828}
@media(max-width:1100px){.feat-grid,.sensor-grid{grid-template-columns:1fr}}
</style>"""


def _build_features_section():
    """Returns HTML for the Features & Logic Overview section."""
    features = [
        ("Schema Auto-Discovery",
         "Automatically identifies time, run/group ID, step, sensor, alarm and SP columns from "
         "CSV headers  --  no hardcoding. Excludes event-counter tallies (arc counts, alarm counts). "
         "Applies a noise filter (within-step std / global range > 13%) and a slow-drift filter "
         "(Spearman > 0.85 in both early/late halves AND very-early range > 20% of global) to "
         "remove sensors that carry no fault-detection signal."),
        ("Inactive-Gas Flow Filter",
         "Each MFC gas token (MFC_AR, MFC_NF3, etc.) is checked: if its primary *Flow column "
         "has median per-run range < 10 sccm, that gas is inactive. ALL sibling columns "
         "(Pressure, Temperature, SP, Ratio) for that gas are cascaded out. A longer active "
         "token (MFC_AR_CLN active) protects its siblings from a shorter inactive token (MFC_AR). "
         "VS_ ratio columns are excluded from driving cascades."),
        ("Aborted Run Detection",
         "Runs with row count < 75% of the median for that tool are classified as Aborted and "
         "excluded from all baseline building, anomaly scoring and step importance analysis. "
         "Shown as orange triangles on trend charts and orange ABORTED badges in the run table."),
        ("Stable-Core 2-Pass Baseline",
         "Pass 1: rough baseline from all runs, score every run. Pass 2: keep bottom 50% "
         "(stable core), rebuild clean baseline excluding the fault bloc. Prevents a persistent "
         "fault (e.g. 20 runs with DCBias at +12 V) from pulling the median toward the fault "
         "value. Early-drift sensors (warmup/consumable) excluded from rough scoring pass only."),
        ("Z-Score Anomaly Scoring & LOO",
         "Score = sum over sensors of exceedance_fraction x mean_excess_z. Threshold = "
         "median + 2.5 x max(MAD, median*0.10, 0.05). Per anomalous run a Leave-One-Out "
         "baseline is rebuilt so the flagged run is never in its own reference. "
         "A continuity check demotes isolated flagged runs unless score > 3x threshold."),
        ("Lead / Lag Divergence Detection",
         "Per anomalous run: find the recipe-elapsed time each sensor first sustains "
         "MIN_CONSEC_DIVERGE (3) consecutive points above 3 sigma using a 2-second rolling "
         "median. Sensors ordered by divergence time; ties broken by mean excess-z. Only the "
         "final recipe step is analysed  --  earlier steps use different setpoints."),
        ("SP Pair Detection & 2-Factor Classification",
         "Auto-detects (actual, SP) column pairs via suffix, prefix, w/r convention, and "
         "infix naming. Each deviating sensor classified as: Genuine Fault (far from baseline "
         "AND far from SP), Change in SP (far from baseline but tracking SP within 5%), or "
         "Anomalous no SP (no setpoint column available)."),
        ("Step Importance Ranking",
         "Per step x sensor: compute per-run mean, then cross-run CV = std / |mean|. Steps "
         "ranked by mean CV across active sensors. Active filter: step mean >= global_min + "
         "5% of global range. Noise floor: std >= 1% of global range. Min 500 rows and "
         "5 runs required. Top-N steps charted with trend charts and correlation heatmaps."),
        ("Highflier Detection & Classification",
         "Runs beyond median +/- 3xIQR per step/sensor are flagged. PATH A: if DATA_QUALITY "
         "varies (IQR > 0), DQ < 50% of median DQ = artifact, else real. PATH B: if DQ is "
         "flat (unreliable), fallback: truncated step data = artifact; consecutive cluster = "
         "real; >= 2 sensors simultaneously elevated = real; isolated spike = artifact."),
        ("Correlation Heatmaps",
         "Per top step: Pearson correlation matrix of per-run means. Aborted runs and "
         "data-artifact highflier runs excluded. Red = strong positive, blue = strong "
         "negative, white = no correlation. Reveals co-moving sensors (e.g. DCBias and "
         "ImpedanceR in an RF match fault) to confirm root cause."),
    ]
    parts = ['<section id="features">',
             '<h2>1 &mdash; Features &amp; Logic Overview</h2>',
             '<p style="font-size:12px;color:#666;margin-bottom:10px">'
             'Click any card to expand the design rationale.</p>',
             '<div class="feat-grid">']
    for title, desc in features:
        parts.append(
            f'<details class="feat-card"><summary>{_esc(title)}</summary>'
            f'<p style="margin-top:8px">{_esc(desc)}</p></details>')
    parts.append('</div></section>')
    return '\n'.join(parts)


def _build_overview_section(tool_results):
    """
    tool_results : { tool_id: { run_id: {rows, score, is_anom, is_abort,
                                          lead_sensor, top_sensors, arc_str} } }
    """
    total   = sum(len(rd) for rd in tool_results.values())
    n_anom  = sum(1 for rd in tool_results.values()
                  for i in rd.values() if i["is_anom"])
    n_abort = sum(1 for rd in tool_results.values()
                  for i in rd.values() if i["is_abort"])
    n_norm  = total - n_anom - n_abort

    parts = ['<section id="overview">',
             '<h2>2 &mdash; All-Runs Status Overview</h2>',
             '<div class="card"><table>'
             f'<tr><th>Total runs</th><td>{total}</td></tr>'
             f'<tr><th>Normal</th><td style="color:#2E7D32"><b>{n_norm}</b></td></tr>'
             f'<tr><th>Anomalous</th><td style="color:#C62828"><b>{n_anom}</b></td></tr>'
             f'<tr><th>Aborted</th><td style="color:#E65100"><b>{n_abort}</b></td></tr>'
             '</table></div>']

    for tool_id, run_info in sorted(tool_results.items()):
        parts.append(f'<h3>Tool {_esc(tool_id)}</h3>')
        parts.append('<div class="card"><table><tr><th>Run ID</th><th>Rows</th>'
                     '<th>Score</th><th>Status</th><th>Lead Sensor</th>'
                     '<th>First ARC</th><th>Top sensors</th></tr>')
        def _run_sort_key(r):
            try: return (0, int(r))
            except (ValueError, TypeError): return (1, r)
        for rid in sorted(run_info, key=_run_sort_key):
            info = run_info[rid]
            if info["is_abort"]:
                chip = '<span class="chip-ab">ABORTED</span>'
                parts.append(
                    f'<tr class="abort"><td>{_esc(rid)}</td>'
                    f'<td>{info["rows"]}</td><td> -- </td><td>{chip}</td>'
                    f'<td></td><td></td>'
                    f'<td style="font-size:11px;color:#E65100">Run ended early</td></tr>')
            elif info["is_anom"]:
                chip = '<span class="chip-a">ANOMALOUS</span>'
                parts.append(
                    f'<tr class="anom"><td>{_esc(rid)}</td>'
                    f'<td>{info["rows"]}</td>'
                    f'<td>{info["score"]:.3f}</td><td>{chip}</td>'
                    f'<td><b>{_esc(info["lead_sensor"])}</b></td>'
                    f'<td>{info["arc_str"]}</td>'
                    f'<td style="font-size:11px;color:#555">{_esc(info["top_sensors"])}</td></tr>')
            else:
                chip = '<span class="chip-n">Normal</span>'
                parts.append(
                    f'<tr class="norm"><td>{_esc(rid)}</td>'
                    f'<td>{info["rows"]}</td>'
                    f'<td>{info["score"]:.3f}</td><td>{chip}</td>'
                    f'<td></td><td></td>'
                    f'<td style="font-size:11px;color:#555">{_esc(info["top_sensors"])}</td></tr>')
        parts.append('</table></div>')
    parts.append('</section>')
    return '\n'.join(parts)


def _build_zscore_section(rep_run_id, rep_tool_id, rep_chain, rep_loo_bl,
                           tool_runs, alarm_info,
                           sensor_classes, sp_pairs, sp_stats, _SENSORS):
    """One z-score chart for the worst anomalous run."""
    if rep_run_id is None:
        return ('<section id="zscore"><h2>3 &mdash; Representative Z-Score</h2>'
                '<div class="card"><p><em>No anomalous runs detected.</em></p>'
                '</div></section>')

    sp_col_names = {sp_col for _, sp_col in (sp_pairs or [])}

    def _is_genuine(sensor, rid):
        if sensor.endswith("_SP") or sensor in sp_col_names:
            return False
        cat = sensor_classes.get(rid, {}).get(sensor, {}).get("category", "")
        if cat in ("Change in SP", "Normal"):
            return False
        if sp_pairs:
            sp_map_local = {a: sc for a, sc in sp_pairs}
            sc = sp_map_local.get(sensor)
            if sc:
                sd = sp_stats.get(rid, {}).get((sensor, sc), {})
                if sd and sd.get("mean_pct", 1.0) < SP_TRACKING_THRESHOLD:
                    return False
        return True

    genuine_chain = [(s, t) for s, t in rep_chain if _is_genuine(s, rep_run_id)]
    if not genuine_chain:
        genuine_chain = rep_chain
    genuine_lead = genuine_chain[0][0] if genuine_chain else "?"
    run_alarms   = (alarm_info or {}).get(rep_run_id, [])
    rrows        = tool_runs[rep_tool_id][rep_run_id]

    svg = svg_zscore_chart(rep_run_id, rrows, rep_loo_bl,
                            genuine_lead, genuine_chain, rep_tool_id,
                            alarms=run_alarms)
    parts = [
        '<section id="zscore">',
        '<h2>3 &mdash; Representative Z-Score Chart</h2>',
        f'<p style="font-size:12px;color:#666">Run <b>{_esc(rep_run_id)}</b> '
        f'(tool {_esc(rep_tool_id)})  --  highest composite anomaly score. '
        f'Lead sensor: <b style="color:#C62828">{_esc(genuine_lead)}</b>. '
        f'±3σ threshold lines shown; vertical markers = first divergence time.</p>',
        f'<div class="card">{svg}</div>',
        '</section>',
    ]
    return '\n'.join(parts)



# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7b — INCREMENTAL HTML WRITER
# Each write_* function receives an open file handle and flushes immediately
# so the browser can display sections as they complete.
# ══════════════════════════════════════════════════════════════════════════════

def _html_header(html_out):
    """Open the HTML file and write the page header + features section."""
    fh = open(html_out, 'w', encoding='utf-8')
    fh.write(
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>"
        f"<title>Combined Analysis Report</title>{HTML_STYLE}</head><body>\n"
        "<h1>Combined Analysis Report &mdash; Recipe Trace + Step Importance</h1>"
        "<p><em>Pure Python &mdash; zero external dependencies.</em></p>\n")
    # Navigation bar
    fh.write(
        '<p style="font-size:12px;background:#E3F2FD;padding:8px 14px;border-radius:4px">'
        + ' &nbsp;| '.join([
            '<a href="#features">1 Features</a>',
            '<a href="#overview">2 Run Overview</a>',
            '<a href="#zscore">3 Z-Score</a>',
            '<a href="#cascade">4 Lead/Lag Cascade</a>',
            '<a href="#envelope">5 Envelope Traces</a>',
            '<a href="#rawtraces">6 Raw Traces</a>',
            '<a href="#spclass">7 SP Classification</a>',
            '<a href="#steps">8 Step Importance</a>',
        ]) + '</p>\n')
    fh.write(_build_features_section() + '\n')
    fh.flush()
    return fh


def _write_anomaly_sections(fh, tool_results, tool_runs, step_info_per_run, looped,
                             alarm_info, sp_pairs, sp_stats, sensor_classes,
                             aborted_all, row_counts_all,
                             lead_lag_all, loo_baselines_all, composite_all_flat,
                             per_sensor_sc_all, anom_set_all, baselines_full_all,
                             rep_run_id, rep_tool_id, sensor_cols):
    """Write sections 2-7 (anomaly detection results) to open file handle."""
    _SENSORS     = sensor_cols or []
    sp_pairs     = sp_pairs or []
    sp_col_names = {sc for _, sc in sp_pairs}

    def _is_genuine(sensor, rid):
        if sensor.endswith("_SP") or sensor in sp_col_names:
            return False
        cat = sensor_classes.get(rid, {}).get(sensor, {}).get("category", "")
        if cat in ("Change in SP", "Normal"):
            return False
        if sp_pairs:
            sp_map_l = {a: sc for a, sc in sp_pairs}
            sc = sp_map_l.get(sensor)
            if sc:
                sd = sp_stats.get(rid, {}).get((sensor, sc), {})
                if sd and sd.get("mean_pct", 1.0) < SP_TRACKING_THRESHOLD:
                    return False
        return True

    # ── Section 2: All-Runs Overview ─────────────────────────────────────────
    print("  Writing section 2 (run overview) ...")
    fh.write(_build_overview_section(tool_results) + '\n')
    fh.flush()

    # ── Section 2b: Primary fault sensor envelope (all runs) ─────────────────
    rep_chain_raw = []
    if rep_run_id:
        rep_chain_raw = lead_lag_all.get(rep_tool_id, {}).get(rep_run_id, (None, [], []))[1]
    genuine_top = [(s, t) for s, t in rep_chain_raw if _is_genuine(s, rep_run_id)] if rep_run_id else []
    primary_sensor = genuine_top[0][0] if genuine_top else (rep_chain_raw[0][0] if rep_chain_raw else None)

    if primary_sensor:
        for tool_id in sorted(tool_runs):
            aborted_set_t  = {rid for tid, rid in aborted_all if tid == tool_id}
            run_dict_clean = {rid: rrows for rid, rrows in tool_runs[tool_id].items()
                              if (tool_id, rid) not in aborted_all}
            run_dict_chart = dict(tool_runs[tool_id])
            anom_set_t  = {rid for t_id, rid in anom_set_all if t_id == tool_id}
            run_ids     = sorted(run_dict_clean)
            run_colors  = {rid: RUN_COLORS[i % len(RUN_COLORS)] for i, rid in enumerate(run_ids)}
            for rid in aborted_set_t:
                run_colors[rid] = "#FF6D00"
            looped_t    = {rid: v for rid, v in looped.items() if rid in run_dict_chart}
            svg_env = svg_envelope_trace(primary_sensor, tool_id, run_dict_chart,
                                          step_info_per_run, looped_t, run_colors,
                                          anom_set_t, alarm_info=alarm_info,
                                          aborted_runs=aborted_set_t)
            if svg_env:
                fh.write(
                    f'<div class="card">'
                    f'<h3>Primary Fault Sensor &mdash; {_esc(primary_sensor)}</h3>'
                    f'<p style="font-size:12px;color:#666">Green band = normal-run envelope '
                    f'(median &plusmn;3&times;MAD). Thick red = anomalous runs. '
                    f'Orange = aborted runs.</p>'
                    f'{svg_env}</div>\n')
        fh.flush()

    # ── Section 3: Representative z-score ────────────────────────────────────
    print("  Writing section 3 (z-score chart) ...")
    if rep_run_id:
        rep_loo_bl    = loo_baselines_all.get(rep_tool_id, {}).get(rep_run_id, {})
        fh.write(_build_zscore_section(
            rep_run_id, rep_tool_id, rep_chain_raw, rep_loo_bl,
            tool_runs, alarm_info, sensor_classes, sp_pairs, sp_stats, _SENSORS) + '\n')
    else:
        fh.write('<section id="zscore"><h2>3 &mdash; Representative Z-Score</h2>'
                 '<div class="card"><p><em>No anomalous runs detected.</em></p>'
                 '</div></section>\n')
    fh.flush()

    # ── Section 4: Lead/lag cascade ──────────────────────────────────────────
    print("  Writing section 4 (lead/lag cascade) ...")
    fh.write('<section id="cascade"><h2>4 &mdash; Lead / Lag Cascade Detail</h2>\n')
    for tool_id in sorted(tool_runs):
        lead_lag = lead_lag_all.get(tool_id, {})
        if not lead_lag:
            continue
        fh.write(f'<h3>Tool {_esc(tool_id)}</h3><div class="card">\n')
        for rid, (lead, chain, non_div) in sorted(lead_lag.items()):
            genuine_chain = [(s, t) for s, t in chain if _is_genuine(s, rid)]
            genuine_lead  = genuine_chain[0][0] if genuine_chain else lead
            lead_t        = genuine_chain[0][1] if genuine_chain else (chain[0][1] if chain else 0)
            run_alarms    = (alarm_info or {}).get(rid, [])
            cascade_names = [s for s, _ in genuine_chain]
            timeline = []
            for ci, (s, t) in enumerate(genuine_chain):
                tag   = "LEAD" if s == genuine_lead else f"lag-{ci}"
                delta = f"+{t - lead_t:.2f}s" if ci > 0 else "0s"
                timeline.append(("sensor", t, tag, s, delta, ""))
            for a in run_alarms:
                delta = f"+{a['elapsed_s'] - lead_t:.2f}s"
                short = a["event_name"].split(".")[-1]
                rel   = classify_alarm_relevance(a, genuine_lead, cascade_names)
                timeline.append(("alarm", a["elapsed_s"], "ALARM", short, delta,
                                  a["description"]))
            timeline.sort(key=lambda x: x[1])
            fh.write(
                f'<p style="margin-top:10px"><strong>Run {_esc(rid)}</strong> &mdash; '
                f'Lead: <span style="color:#C62828;font-weight:bold">'
                f'{_esc(genuine_lead)}</span></p>\n')
            fh.write('<table><tr><th>Time (s)</th><th>Type</th>'
                     '<th>Sensor/Event</th><th>Delta</th><th>Detail</th></tr>\n')
            for entry in timeline:
                kind = entry[0]
                if kind == "sensor":
                    _, t, tag, s, delta, _ = entry
                    col = "#C62828" if tag == "LEAD" else "#E65100"
                    fh.write(f'<tr><td>{t:.1f}</td>'
                             f'<td><b style="color:{col}">{_esc(tag)}</b></td>'
                             f'<td>{_esc(s)}</td><td>{_esc(delta)}</td><td></td></tr>\n')
                else:
                    _, t, tag, short, delta, desc = entry[:6]
                    fh.write(f'<tr><td>{t:.1f}</td><td>{_esc(tag)}</td>'
                             f'<td>{_esc(short)}</td><td>{_esc(delta)}</td>'
                             f'<td style="font-size:11px;color:#555">'
                             f'{_esc(str(desc)[:80])}</td></tr>\n')
            fh.write('</table>\n')
        fh.write('</div>\n')
    fh.write('</section>\n')
    fh.flush()

    # ── Section 5: Baseline envelope traces ──────────────────────────────────
    print("  Writing section 5 (envelope traces) ...")
    fh.write('<section id="envelope"><h2>5 &mdash; Baseline Envelope Traces</h2>\n')
    for tool_id in sorted(tool_runs):
        lead_lag       = lead_lag_all.get(tool_id, {})
        anom_set_t     = {rid for t_id, rid in anom_set_all if t_id == tool_id}
        aborted_set_t  = {rid for tid, rid in aborted_all if tid == tool_id}
        run_dict_clean = {rid: rrows for rid, rrows in tool_runs[tool_id].items()
                          if (tool_id, rid) not in aborted_all}
        # Chart dict includes aborted runs so they appear as orange traces
        run_dict_chart = dict(tool_runs[tool_id])
        run_ids        = sorted(run_dict_clean)
        run_colors     = {rid: RUN_COLORS[i % len(RUN_COLORS)] for i, rid in enumerate(run_ids)}
        for rid in aborted_set_t:
            run_colors[rid] = "#FF6D00"
        looped_t       = {rid: v for rid, v in looped.items() if rid in run_dict_chart}
        per_sensor_sc  = per_sensor_sc_all.get(tool_id, {})

        env_sensors = []; seen_env = set()
        for rid, (lead, chain, _) in sorted(lead_lag.items()):
            for s, _ in chain:
                if _is_genuine(s, rid) and s not in seen_env:
                    env_sensors.append(s); seen_env.add(s)
        for rid in sorted(anom_set_t):
            for s, sc_v in sorted(per_sensor_sc.get(rid, {}).items(), key=lambda x: -x[1]):
                if s not in seen_env and sc_v > 0.0 and _is_genuine(s, rid):
                    env_sensors.append(s); seen_env.add(s)

        if env_sensors:
            fh.write(f'<h3>Tool {_esc(tool_id)}</h3><div class="card">\n')
            for sensor in env_sensors:
                has = any(_f(r.get(sensor)) is not None
                          for rrows in run_dict_chart.values() for r in rrows)
                if not has:
                    continue
                svg = svg_envelope_trace(sensor, tool_id, run_dict_chart,
                                          step_info_per_run, looped_t, run_colors,
                                          anom_set_t, alarm_info=alarm_info,
                                          aborted_runs=aborted_set_t)
                if svg:
                    fh.write(svg + '\n')
            fh.write('</div>\n')
    fh.write('</section>\n')
    fh.flush()

    # ── Section 6: Raw sensor traces ─────────────────────────────────────────
    print("  Writing section 6 (raw traces) ...")
    fh.write('<section id="rawtraces"><h2>6 &mdash; Raw Sensor Traces</h2>\n')
    for tool_id in sorted(tool_runs):
        lead_lag       = lead_lag_all.get(tool_id, {})
        anom_set_t     = {rid for t_id, rid in anom_set_all if t_id == tool_id}
        aborted_set_t  = {rid for tid, rid in aborted_all if tid == tool_id}
        run_dict_clean = {rid: rrows for rid, rrows in tool_runs[tool_id].items()
                          if (tool_id, rid) not in aborted_all}
        # Chart dict includes aborted runs so they appear as orange traces
        run_dict_chart = dict(tool_runs[tool_id])
        run_ids        = sorted(run_dict_clean)
        run_colors     = {rid: RUN_COLORS[i % len(RUN_COLORS)] for i, rid in enumerate(run_ids)}
        for rid in aborted_set_t:
            run_colors[rid] = "#FF6D00"
        looped_t       = {rid: v for rid, v in looped.items() if rid in run_dict_chart}

        ordered = []; seen_s = set()
        for rid, (lead, chain, _) in sorted(lead_lag.items()):
            for s, _ in chain:
                if _is_genuine(s, rid) and s not in seen_s:
                    ordered.append(s); seen_s.add(s)
        if ordered:
            fh.write(f'<h3>Tool {_esc(tool_id)}</h3><div class="card">\n')
            for sensor in ordered:
                has = any(_f(r.get(sensor)) is not None
                          for rrows in run_dict_chart.values() for r in rrows)
                if not has:
                    continue
                svg = svg_raw_trace(sensor, tool_id, run_dict_chart,
                                     step_info_per_run, looped_t, run_colors,
                                     anom_set_t, alarm_info=alarm_info,
                                     aborted_runs=aborted_set_t)
                if svg:
                    fh.write(svg + '\n')
            fh.write('</div>\n')
    fh.write('</section>\n')
    fh.flush()

    # ── Section 7: SP classification ─────────────────────────────────────────
    print("  Writing section 7 (SP classification) ...")
    fh.write('<section id="spclass"><h2>7 &mdash; Sensor Deviation Classification</h2>\n')
    if sensor_classes:
        CAT_COLOR = {
            "Genuine Fault":     ("#B71C1C", "#FFEBEE"),
            "Anomalous (no SP)": ("#E65100", "#FFF3E0"),
            "Change in SP":      ("#1565C0", "#E3F2FD"),
            "Normal":            ("#2E7D32", "#E8F5E9"),
        }
        for tool_id in sorted(tool_runs):
            anom_set_t = {rid for t_id, rid in anom_set_all if t_id == tool_id}
            if not anom_set_t:
                continue
            fh.write(f'<h3>Tool {_esc(tool_id)}</h3><div class="card">\n')
            fh.write('<table><tr><th>Run</th><th>Sensor</th><th>Category</th>'
                     '<th>Baseline z</th><th>SP deviation</th></tr>\n')
            for rid in sorted(anom_set_t):
                for s, info in sorted(sensor_classes.get(rid, {}).items()):
                    cat = info.get("category", "")
                    if cat in ("Normal", ""):
                        continue
                    tc, bc = CAT_COLOR.get(cat, ("#333", "white"))
                    bl_z   = info.get("bl_z", float("nan"))
                    sp_pct = info.get("sp_pct_dev", None)
                    sp_str = f"{sp_pct*100:.1f}%" if sp_pct is not None else "&mdash;"
                    fh.write(
                        f'<tr style="background:{bc}"><td>{_esc(rid)}</td>'
                        f'<td>{_esc(s)}</td>'
                        f'<td style="color:{tc};font-weight:bold">{_esc(cat)}</td>'
                        f'<td>{bl_z:.2f}</td><td>{sp_str}</td></tr>\n')
            fh.write('</table></div>\n')
            # SP deviation charts for anomalous runs in this tool
            if sp_pairs:
                run_dict_t   = {rid: rrows for rid, rrows in tool_runs[tool_id].items()
                                if (tool_id, rid) not in aborted_all}
                run_ids_t    = sorted(run_dict_t)
                run_colors_t = {rid: RUN_COLORS[i % len(RUN_COLORS)]
                                for i, rid in enumerate(run_ids_t)}
                for actual_col, sp_col in sp_pairs:
                    has_anom = any(
                        sp_stats.get(rid, {}).get((actual_col, sp_col)) is not None
                        for rid in anom_set_t)
                    if not has_anom:
                        continue
                    svg = svg_sp_deviation(actual_col, sp_col, run_dict_t,
                                           sp_stats, run_colors_t, anom_set_t)
                    if svg:
                        fh.write(f'<div class="card">'
                                 f'<h4 style="margin:6px 0">{_esc(actual_col)}'
                                 f' vs Setpoint</h4>{svg}</div>\n')
    else:
        fh.write('<div class="card"><p>No SP pairs detected.</p></div>\n')
    fh.write('</section>\n')
    fh.flush()
    print("  Anomaly sections written and flushed.")


def _write_step_sections(fh, step_info_si, run_order_si, run_tool_si,
                          run_means_by_step, highfliers,
                          aborted_all, row_counts_all,
                          forced_step_sensors=None,
                          noisy_sensor_steps=None):
    """Write section 8 (step importance) to open file handle, then close it."""
    forced_step_sensors  = forced_step_sensors  or {}
    noisy_sensor_steps   = noisy_sensor_steps   or set()
    tools_si  = sorted(set(run_tool_si.values()))
    # Top steps by CV ranking
    cv_top_steps = sorted(step_info_si, key=lambda s: -step_info_si[s]["mean_cv"])[:TOP_N_STEPS]
    # Forced steps: steps not already in top_steps but containing anomaly-flagged sensors
    forced_only_steps = [sn for sn in forced_step_sensors if sn not in cv_top_steps
                         and sn in step_info_si]
    top_steps = cv_top_steps + forced_only_steps
    all_steps = sorted(step_info_si, key=lambda s: -step_info_si[s]["mean_cv"])

    print("  Writing section 8 (step importance ranking) ...")
    fh.write('<section id="steps"><h2>8 &mdash; Step Importance Ranking</h2>\n')
    fh.write(
        f'<p style="font-size:12px;color:#666">Steps ranked by mean cross-run CV '
        f'(std/|mean|). Aborted runs excluded. '
        f'Top {TOP_N_STEPS} steps charted below.</p>\n')

    fh.write('<div class="card"><table><tr><th>Rank</th><th>Step</th><th>Name</th>'
             '<th>Rows</th><th>Runs</th><th>Mean CV</th><th>Top sensors</th></tr>\n')
    for rank, sn in enumerate(all_steps, 1):
        si     = step_info_si[sn]
        is_top    = sn in cv_top_steps
        is_forced = sn in forced_only_steps
        cls    = ' class="top1"' if rank == 1 else ""
        if is_top:
            badge = '<span class="badge-top">CHARTED</span>'
        elif is_forced:
            badge = ('<span class="badge-top" style="background:#D32F2F">'
                     'ANOMALY FLAGGED</span>')
        else:
            badge = ""
        top3   = si["sensor_cvs"][:3]
        top_str = "&nbsp;&nbsp;".join(
            f'<b>{_esc(_short(s))}</b> CV={cv:.3f}' for s, cv, *_ in top3)
        fh.write(
            f'<tr{cls}><td>{rank}{badge}</td><td>Step {_esc(sn)}</td>'
            f'<td>{_esc(si["name"])}</td><td>{si["rows"]:,}</td>'
            f'<td>{si["n_runs"]}</td><td><b>{si["mean_cv"]:.4f}</b></td>'
            f'<td style="font-size:11px">{top_str}</td></tr>\n')
    fh.write('</table></div>\n')
    fh.flush()

    # Per-step detail
    for sn in top_steps:
        si      = step_info_si.get(sn, {})
        if not si:
            continue
        rm_step = run_means_by_step.get(sn, {})
        hf_step = highfliers.get(sn, {})

        is_forced = sn in forced_only_steps
        cv_rank   = cv_top_steps.index(sn) + 1 if sn in cv_top_steps else None
        if is_forced:
            badge = '<span class="badge-top" style="background:#D32F2F">ANOMALY FLAGGED</span>'
            reason_note = ('Step included because anomaly-flagged sensors were detected here '
                           'during per-sensor independent analysis.')
        else:
            badge = f'<span class="badge-top">Top {cv_rank}</span>'
            reason_note = None

        print(f"  Writing step {sn} ({si['name']}) detail ...")
        fh.write(
            f'<h3>Step {_esc(sn)} &mdash; {_esc(si["name"])}'
            f'{badge}</h3>\n'
            f'<p style="font-size:12px;color:#666">{si["rows"]:,} rows, '
            f'{si["n_runs"]} runs. Mean CV = <b>{si["mean_cv"]:.4f}</b>.</p>\n')
        if reason_note:
            fh.write(f'<p style="font-size:12px;color:#D32F2F">'
                     f'&#9679; {_esc(reason_note)}</p>\n')

        # Forced sensors for this step (not already in sensor_cvs)
        cv_sensor_names = {s for s, *_ in si.get("sensor_cvs", [])}
        forced_sensors_this_step = forced_step_sensors.get(sn, set()) - cv_sensor_names

        # Stats table -- one sub-table per tool
        for _t_stat in tools_si:
            _t_run_set = {ck for ck in run_tool_si if run_tool_si[ck] == _t_stat}
            _short_t   = "_".join(_t_stat.split("_")[-3:]) if "_" in _t_stat else _t_stat
            fh.write(f'<div class="card"><h3>Sensor Statistics &mdash; '
                     f'{_esc(_short_t)}</h3>'
                     f'<table><tr><th>Sensor</th><th>Mean</th><th>Std</th>'
                     f'<th>Min run-mean</th><th>Max run-mean</th>'
                     f'<th>Range</th><th>CV</th><th>Highfliers</th></tr>\n')
            max_cv_t = max((cv for s, cv, *_ in si.get("sensor_cvs", [])), default=1e-9)
            # All sensors to show: CV-ranked + forced sensors for this step
            _cv_rows   = [(s, cv, False) for s, cv, *_ in si.get("sensor_cvs", [])]
            _forced_rows = [(s, 0.0, True) for s in sorted(forced_sensors_this_step)]
            for s, cv, _is_forced in _cv_rows + _forced_rows:
                cv = cv  # already set above, re-assign for clarity
                rvals = [v for ck, v in rm_step.get(s, {}).items()
                         if ck in _t_run_set]
                if not rvals:
                    continue
                import math as _math
                rmin  = min(rvals); rmax = max(rvals)
                tmn   = sum(rvals) / len(rvals)
                tstd  = _math.sqrt(sum((v-tmn)**2 for v in rvals)/len(rvals)) if len(rvals)>1 else 0.0
                trng  = rmax - rmin
                tcv   = tstd / abs(tmn) if abs(tmn) > 1e-9 else 0.0
                bar_w = min(int(tcv / max(max_cv_t, 1e-9) * 60), 60)
                bar   = (f'<span style="display:inline-block;width:{bar_w}px;height:8px;'
                         f'background:#1565C0;border-radius:2px;vertical-align:middle;'
                         f'margin-right:4px"></span>')
                hf_s  = {ck: info for ck, info in hf_step.get(s, {}).items()
                         if ck in _t_run_set}
                if hf_s:
                    real_r = [r for r, h in hf_s.items() if h["class"] == "real"]
                    art_r  = [r for r, h in hf_s.items() if h["class"] == "artifact"]
                    hf_str = ""
                    if real_r:
                        shown = ", ".join(sorted(real_r)[:3]) + ("\u2026" if len(real_r) > 3 else "")
                        hf_str += (f'<span style="color:{HIGHFLIER_REAL};font-weight:bold">'
                                   f'{len(real_r)} real</span> '
                                   f'<span style="font-size:10px;color:#888">({shown})</span> ')
                    if art_r:
                        shown = ", ".join(sorted(art_r)[:3]) + ("\u2026" if len(art_r) > 3 else "")
                        hf_str += (f'<span style="color:#9E6C00;font-weight:bold">'
                                   f'{len(art_r)} artifact</span> '
                                   f'<span style="font-size:10px;color:#888">({shown})</span>')
                else:
                    hf_str = '<span style="color:#aaa">none</span>'
                forced_badge = (' <span style="background:#D32F2F;color:white;'
                                'font-size:9px;padding:1px 4px;border-radius:2px">'
                                'ANOMALY</span>') if _is_forced else ''
                fh.write(
                    f'<tr><td><b>{_esc(s)}</b>{forced_badge}<br>'
                    f'<span style="font-size:10px;color:#888">{_esc(_short(s))}</span></td>'
                    f'<td>{tmn:.4g}</td><td>{tstd:.4g}</td>'
                    f'<td>{rmin:.4g}</td><td>{rmax:.4g}</td>'
                    f'<td>{trng:.4g}</td><td>{bar}{tcv:.4f}</td>'
                    f'<td style="font-size:11px">{hf_str}</td></tr>\n')
            fh.write('</table></div>\n')
        fh.flush()

        # Highflier detail
        for _t_hf in tools_si:
            _t_hf_run_set = {ck for ck in run_tool_si if run_tool_si[ck] == _t_hf}
            _short_t_hf   = "_".join(_t_hf.split("_")[-3:]) if "_" in _t_hf else _t_hf
            all_hf = [(s, ckey, info) for s, flags in hf_step.items()
                      for ckey, info in flags.items() if ckey in _t_hf_run_set]
            if not all_hf:
                continue
            fh.write(f'<div class="card"><h3>Highflier Detail &mdash; '
                     f'{_esc(_short_t_hf)}</h3>'
                     f'<p style="font-size:12px;color:#666">Runs beyond '
                     f'median &plusmn;{HIGHFLIER_IQR_MULT}&times;IQR.</p>'
                     '<table><tr><th>Sensor</th><th>Run</th><th>Value</th>'
                     '<th>Dir</th><th>Class</th><th>Reason</th></tr>\n')
            for s, ckey, info in sorted(all_hf, key=lambda x: (x[2]["class"], x[0], x[1])):
                rid_lbl = ckey.split("::", 1)[1] if "::" in ckey else ckey
                cls    = info["class"]
                hf_col = HIGHFLIER_REAL if cls == "real" else HIGHFLIER_ART
                lbl    = "REAL ISSUE" if cls == "real" else "DATA ARTIFACT"
                fh.write(
                    f'<tr class="{cls}"><td>{_esc(s)}</td><td><b>{_esc(rid_lbl)}</b></td>'
                    f'<td>{info["value"]:.4g}</td><td>{info["direction"]}</td>'
                    f'<td><span style="background:{hf_col};color:white;padding:1px 6px;'
                    f'border-radius:3px;font-size:11px">{lbl}</span></td>'
                    f'<td style="font-size:11px">{_esc(info["reason"])}</td></tr>\n')
            fh.write('</table></div>\n')
        fh.flush()

        # Trend charts -- one sub-section per tool so each tool's runs are
        # plotted independently with its own mean+/-std band.
        fh.write('<div class="card"><h3>Per-Run Trend Charts</h3>\n')
        for _tool in tools_si:
            # Filter run_order and run_means to this tool only
            _tool_run_order = [ck for ck in run_order_si
                               if run_tool_si.get(ck) == _tool]
            if not _tool_run_order:
                continue
            _tool_run_set   = set(_tool_run_order)
            _tool_run_tool  = {ck: _tool for ck in _tool_run_order}
            _short_tool     = "_".join(_tool.split("_")[-3:]) if "_" in _tool else _tool

            fh.write(f'<h4 style="margin:10px 0 4px;color:#37474F">'
                     f'Tool: {_esc(_short_tool)}</h4>'
                     f'<div class="sensor-grid">\n')

            # CV-ranked sensors
            for s, cv, mn, std, rng in si.get("sensor_cvs", []):
                rm_sensor_all  = rm_step.get(s, {})
                rm_sensor_tool = {ck: v for ck, v in rm_sensor_all.items()
                                  if ck in _tool_run_set}
                if len(rm_sensor_tool) < 3:
                    continue
                hf_sensor = {ck: info for ck, info in hf_step.get(s, {}).items()
                             if ck in _tool_run_set}
                svg = svg_trend_chart(s, sn, rm_sensor_tool,
                                      _tool_run_order, _tool_run_tool,
                                      [_tool],
                                      aborted_runs=aborted_all,
                                      row_counts=row_counts_all,
                                      highflier_info=hf_sensor)
                if svg:
                    fh.write(f'<div>{svg}</div>\n')

            # Forced sensors (anomaly-flagged, not in CV-ranked list)
            if forced_sensors_this_step:
                fh.write('</div>\n'
                         '<p style="font-size:11px;color:#D32F2F;margin:4px 0">'
                         '&#9679; Sensors below were not top-ranked by CV but had '
                         'anomalous runs detected by per-sensor analysis:</p>'
                         '<div class="sensor-grid">\n')
                for s in sorted(forced_sensors_this_step):
                    rm_sensor_all  = rm_step.get(s, {})
                    rm_sensor_tool = {ck: v for ck, v in rm_sensor_all.items()
                                      if ck in _tool_run_set}
                    if len(rm_sensor_tool) < 3:
                        continue
                    svg = svg_trend_chart(s, sn, rm_sensor_tool,
                                          _tool_run_order, _tool_run_tool,
                                          [_tool],
                                          aborted_runs=aborted_all,
                                          row_counts=row_counts_all)
                    if svg:
                        fh.write(f'<div>{svg}</div>\n')

            fh.write('</div>\n')
        fh.write('</div>\n')
        fh.flush()

        # Correlation heatmap -- one per tool so HW differences do not mix
        artifact_runs = {rid for flags in hf_step.values()
                         for rid, info in flags.items() if info.get("class") == "artifact"}
        excl_corr = aborted_all | artifact_runs
        # Exclude sensors that are flagged as noisy for this specific step --
        # their run-means are noise-driven and would pollute the correlations.
        rm_step_clean = {s: rd for s, rd in rm_step.items()
                         if (s, sn) not in noisy_sensor_steps}
        # collect all composite keys present in rm_step_clean, grouped by tool
        _all_rm_runs = set()
        for _sv in rm_step_clean.values():
            _all_rm_runs.update(_sv.keys())
        _tool_run_ids = defaultdict(list)
        for _ckey in _all_rm_runs:
            _tool_run_ids[run_tool_si.get(_ckey, "ALL")].append(_ckey)
        # Always render one heatmap per tool (never skip a tool)
        _any_corr = False
        for _t_name in sorted(_tool_run_ids):
            # Exclude runs from all other tools so each heatmap is tool-specific
            _excl_this = excl_corr | (_all_rm_runs - set(_tool_run_ids[_t_name]))
            _cs, _cm = compute_correlations(rm_step_clean, _excl_this)
            if not _cs:
                continue
            if not _any_corr:
                fh.write('<div class="card"><h3>Sensor Correlation Matrix</h3>\n')
                _any_corr = True
            fh.write(f'<h4 style="margin:12px 0 4px">'
                     f'Tool: {_esc(_t_name)}</h4>\n')
            n_runs_t = len(_tool_run_ids[_t_name]) - len(excl_corr & set(_tool_run_ids[_t_name]))
            fh.write(f'<p style="font-size:12px;color:#666">Pearson r of per-run means '
                     f'({n_runs_t} runs). '
                     f'<span style="color:#C62828">Red</span>=positive, '
                     f'<span style="color:#1565C0">Blue</span>=negative.</p>\n')
            fh.write(svg_corr_heatmap(_cs, _cm, si["name"]) + '\n')
        if _any_corr:
            fh.write('</div>\n')
            fh.flush()

    fh.write('</section>\n')

    # Footer + close
    fh.write(
        '<div style="text-align:center;font-size:11px;color:#aaa;margin-top:40px">'
        'recipe_analysis_report.py &mdash; zero external dependencies'
        '</div></body></html>\n')
    fh.flush()
    fh.close()
    print("  Step sections written. File closed.")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8 — MAIN
# ══════════════════════════════════════════════════════════════════════════════

def _prompt_data_file():
    print("\n" + "=" * 70)
    print("  Combined Analysis Report")
    print("=" * 70)
    while True:
        raw = input("\n  Enter path to CSV file or directory:\n  > ").strip().strip('"').strip("'")
        if os.path.isfile(raw):
            data_path = raw; break
        if os.path.isdir(raw):
            files = sorted(f for f in os.listdir(raw)
                           if os.path.splitext(f)[1].lower() in (".csv", ".xlsx", ".xlsm"))
            if not files:
                print("  No CSV files found."); continue
            for i, f in enumerate(files):
                print(f"  {i+1}. {f}")
            while True:
                try:
                    idx = int(input("  Pick number: ").strip()) - 1
                    data_path = os.path.join(raw, files[idx]); break
                except (ValueError, IndexError):
                    print("  Invalid.")
            break
        print(f"  Not found: {raw}")
    stem     = os.path.splitext(data_path)[0]
    html_out = stem + "_combined_report.html"
    return data_path, html_out


def main():
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    filepath, html_out = _prompt_data_file()

    # ── 1. Load ───────────────────────────────────────────────────────────────
    print(f"\nLoading {filepath} ...")
    rows = load_csv(filepath)
    print(f"  {len(rows):,} rows loaded.")
    headers = list(rows[0].keys()) if rows else []

    # ── 2. Schema discovery ───────────────────────────────────────────────────
    print("\nDiscovering schema ...")
    schema      = discover_schema(rows, headers)
    sensor_cols = schema.get("sensor_cols", [])
    if not sensor_cols:
        print("  ERROR: no numeric sensor columns detected."); sys.exit(1)
    print(f"  Sensors auto-detected: {len(sensor_cols)}")

    _flush_stdin()
    raw_ex = input(
        "\n  Additional sensor columns to EXCLUDE (space/comma separated)\n"
        "  or Enter to skip:\n  > ").strip()
    EXCLUDED_SENSORS.clear()
    if raw_ex:
        for tok in re.split(r'[\s,;]+', raw_ex):
            if tok:
                EXCLUDED_SENSORS.add(tok.strip())
    sensor_cols = [s for s in sensor_cols
                   if not any(ex.lower() in s.lower() or s.lower() in ex.lower()
                               for ex in EXCLUDED_SENSORS)]
    print(f"  Sensors after exclusions: {len(sensor_cols)}")

    _flush_stdin()
    raw_ka = input(
        "\n  KNOWN ANOMALOUS run IDs (space/comma) or Enter for auto-detection:\n"
        "  > ").strip()
    KNOWN_ANOMALOUS.clear()
    if raw_ka:
        for tok in re.split(r'[\s,;]+', raw_ka):
            if tok:
                KNOWN_ANOMALOUS.add(tok.strip())

    rough_scoring_excl = schema.get("rough_scoring_excl", set())

    # ── Open HTML file and write header + features immediately ────────────────
    print(f"\nOpening output file: {html_out}")
    print("  (You can open this file in a browser now — sections are written as they complete)")
    fh = _html_header(html_out)
    print("  Features section written.")

    # ── 3. Prepare runs ───────────────────────────────────────────────────────
    print("\nPreparing runs ...")
    tool_runs, step_info_per_run, alarm_info = prepare_runs(rows, schema)
    if not tool_runs:
        schema2 = dict(schema); schema2["step_num_col"] = None
        tool_runs, step_info_per_run, alarm_info = prepare_runs(rows, schema2)
    if not tool_runs:
        print("  ERROR: no run data found."); sys.exit(1)
    looped = detect_looped_steps(step_info_per_run)

    # ── 4. Aborted run detection ──────────────────────────────────────────────
    print("\nDetecting aborted runs ...")
    aborted_all, row_counts_all = detect_aborted_runs(tool_runs)
    if not aborted_all:
        print("  None detected.")

    # ── 5. SP pairs ───────────────────────────────────────────────────────────
    sp_pairs = detect_sp_pairs(rows)
    print(f"  SP pairs detected: {len(sp_pairs)}")

    # ── 6. Recipe-trace anomaly analysis ─────────────────────────────────────
    print("\nRunning anomaly analysis (phase 1 of 2) ...")
    _SENSORS = sensor_cols

    tool_results       = {}
    lead_lag_all       = {}
    loo_baselines_all  = {}
    composite_all_flat = {}
    per_sensor_sc_all  = {}
    anom_set_all       = set()
    baselines_full_all = {}
    sensor_classes_all = {}
    sp_stats_all       = {}
    rep_run_id = rep_tool_id = None
    best_score = -1.0

    for tool_id in sorted(tool_runs):
        run_dict       = tool_runs[tool_id]
        run_dict_clean = {rid: rrows for rid, rrows in run_dict.items()
                          if (tool_id, rid) not in aborted_all}
        run_ids        = sorted(run_dict)
        use_stable_core = len(run_dict_clean) >= 10

        baselines_rough = {}
        for s in _SENSORS:
            bl = build_baseline(run_dict_clean, s)
            if bl[0]:
                baselines_rough[s] = bl

        if use_stable_core:
            rough_scores = {}
            for rid, rrows in run_dict_clean.items():
                rough_scores[rid] = sum(
                    compute_run_score(rrows, bl, s)
                    for s, bl in baselines_rough.items()
                    if s not in rough_scoring_excl)
            score_vals  = sorted(rough_scores.values())
            score_med   = score_vals[len(score_vals) // 2]
            stable_core = {rid for rid, sc in rough_scores.items() if sc <= score_med}
            if len(stable_core) < 3:
                stable_core = set(run_dict_clean.keys())
            baselines_full = {}
            for s in _SENSORS:
                bl = build_baseline(run_dict_clean, s,
                                    exclude=set(run_dict_clean) - stable_core)
                if bl[0]:
                    baselines_full[s] = bl
            for s, bl in baselines_rough.items():
                if s not in baselines_full:
                    baselines_full[s] = bl
        else:
            baselines_full = baselines_rough
            stable_core    = set(run_dict_clean.keys())

        baselines_full_all[tool_id] = baselines_full

        # Precompute MFC ramp intervals for this tool (Feature 3)
        mfc_ramp_ivs = {}
        for s in _SENSORS:
            if _is_mfc_flow_sensor(s):
                ivs = _build_mfc_ramp_intervals(run_dict_clean, s)
                if ivs:
                    mfc_ramp_ivs[s] = ivs

        composite = {}; per_sensor_sc = defaultdict(dict)
        for rid, rrows in run_dict_clean.items():
            total = 0.0
            for s, bl in baselines_full.items():
                sc = compute_run_score(rrows, bl, s,
                                       ramp_intervals=mfc_ramp_ivs.get(s))
                per_sensor_sc[rid][s] = sc; total += sc
            composite[rid] = total

        # Per-step scoring pass: score each step independently so a
        # localised shift confined to one step (e.g. Dep DCBias) is
        # not diluted by noise across the other steps.
        step_nums_in_data = set()
        for _rrows_v in run_dict_clean.values():
            for _r in _rrows_v:
                _sn_v = _r.get("_step_number", "")
                if _sn_v: step_nums_in_data.add(_sn_v)
        step_nums_in_data -= SKIP_STEPS
        for _sn_v in sorted(step_nums_in_data):
            _sf = {_sn_v}
            _step_comp = {}
            for rid, rrows in run_dict_clean.items():
                _step_comp[rid] = sum(
                    compute_run_score(rrows, bl, s, step_filter=_sf,
                                      ramp_intervals=mfc_ramp_ivs.get(s))
                    for s, bl in baselines_full.items())
            _step_anom = flag_anomalous(_step_comp)
            for rid, _is_a in _step_anom.items():
                if _is_a:
                    composite[rid] = max(composite.get(rid, 0.0),
                                         _step_comp[rid])

        # Piecewise anomaly detection (Feature 5): if one sensor dominates
        # the composite score, segregate it and score independently to
        # ensure runs with smaller-but-real anomalies are still detected.
        segregated_sensors = _segregate_outlier_sensors(per_sensor_sc)
        if segregated_sensors:
            print(f"  Tool {tool_id}: segregating outlier sensor(s) for independent "
                  f"scoring: {sorted(segregated_sensors)}")
            composite_normal = {}
            composite_segregated = {}
            for rid in run_dict_clean:
                total_normal = 0.0
                total_seg    = 0.0
                for s, bl in baselines_full.items():
                    sc = per_sensor_sc[rid].get(s, 0.0)
                    if s in segregated_sensors:
                        total_seg    += sc
                    else:
                        total_normal += sc
                composite_normal[rid]     = total_normal
                composite_segregated[rid] = total_seg
            anomalous_normal = flag_anomalous(composite_normal) if composite_normal else {}
            anomalous_seg    = flag_anomalous(composite_segregated) if composite_segregated else {}
            # Merge: boost composite score for runs flagged in either pool
            for rid in composite:
                if anomalous_seg.get(rid, False):
                    composite[rid] = max(composite[rid],
                                         composite_segregated.get(rid, 0.0))

        anomalous  = flag_anomalous(composite)
        anom_set_t = {rid for rid, a in anomalous.items() if a}
        anom_set_all |= {(tool_id, rid) for rid in anom_set_t}
        composite_all_flat.update(composite)
        per_sensor_sc_all[tool_id] = dict(per_sensor_sc)

        lead_lag = {}; loo_baselines = {}
        for rid in sorted(anom_set_t):
            if use_stable_core and rid not in stable_core:
                loo_bl = baselines_full
            else:
                excl = (set(run_dict_clean) - stable_core) | {rid}
                loo_bl = {}
                for s in _SENSORS:
                    bl = build_baseline(run_dict_clean, s, exclude=excl)
                    if bl[0]:
                        loo_bl[s] = bl
                for s, bl in baselines_full.items():
                    if s not in loo_bl:
                        loo_bl[s] = bl
            loo_baselines[rid] = loo_bl
            lead, chain, non_div = detect_lead_lag(run_dict_clean[rid], loo_bl, _SENSORS,
                                                   mfc_ramp_ivs=mfc_ramp_ivs)
            if lead:
                lead_lag[rid] = (lead, chain, non_div)
            if composite.get(rid, 0) > best_score:
                best_score = composite[rid]; rep_run_id = rid; rep_tool_id = tool_id

        lead_lag_all[tool_id]      = lead_lag
        loo_baselines_all[tool_id] = loo_baselines

        sc_t = sp_st = {}
        if sp_pairs:
            sc_t, sp_st = classify_sensor_deviations(
                run_dict_clean, sp_pairs, anom_set_t, sensor_list=_SENSORS)
        sensor_classes_all.update(sc_t); sp_stats_all.update(sp_st)

        # build tool_results for overview
        sp_col_names_t = {sp_col for _, sp_col in sp_pairs}
        tool_results[tool_id] = {}
        for rid in run_ids:
            is_abort = (tool_id, rid) in aborted_all
            is_anom  = rid in anom_set_t
            sc_val   = composite.get(rid, 0.0)
            lead_str = ""
            if is_anom and rid in lead_lag:
                lchain = lead_lag[rid][1]
                def _is_gen(s):
                    if s.endswith("_SP") or s in sp_col_names_t: return False
                    cat = sc_t.get(rid, {}).get(s, {}).get("category", "")
                    return cat not in ("Change in SP", "Normal")
                gc = [(s, t) for s, t in lchain if _is_gen(s)]
                lead_str = gc[0][0] if gc else (lchain[0][0] if lchain else "")
            top_s   = sorted(per_sensor_sc.get(rid, {}).items(), key=lambda x: -x[1])
            top_str = ", ".join(f"{s}({v:.2f})" for s, v in top_s[:3] if v > 0.001)
            arc_ev  = [a for a in (alarm_info or {}).get(rid, []) if a.get("is_arc_onset")]
            if arc_ev:
                fa = min(arc_ev, key=lambda a: a["elapsed_s"])
                arc_str = (f'<span style="color:#B71C1C;font-weight:bold">'
                           f'T={fa["elapsed_s"]:.1f}s</span>'
                           f'<br><span style="font-size:10px;color:#777">'
                           f'{fa.get("arc_col","").split("_")[-1]}</span>')
            else:
                arc_str = '<span style="color:#aaa">&mdash;</span>'
            tool_results[tool_id][rid] = dict(
                rows=len(run_dict[rid]), score=sc_val,
                is_anom=is_anom, is_abort=is_abort,
                lead_sensor=lead_str, top_sensors=top_str, arc_str=arc_str)

        print(f"  Tool {tool_id}: {len(run_dict_clean)} runs, "
              f"{len(anom_set_t)} anomalous")

    # ── Write anomaly sections to HTML immediately ────────────────────────────
    print("\nWriting anomaly sections to HTML (phase 1 complete) ...")
    _write_anomaly_sections(
        fh, tool_results, tool_runs, step_info_per_run, looped,
        alarm_info, sp_pairs, sp_stats_all, sensor_classes_all,
        aborted_all, row_counts_all,
        lead_lag_all, loo_baselines_all, composite_all_flat,
        per_sensor_sc_all, anom_set_all, baselines_full_all,
        rep_run_id, rep_tool_id, sensor_cols)
    print(f"  -> Sections 2-7 visible in browser: {html_out}")

    # ── 7. Step importance analysis (phase 2) ─────────────────────────────────
    print("\nRunning step importance analysis (phase 2 of 2) ...")
    run_signals, _, median_step_rows = build_run_process_signals(
        rows, schema["run_col"], schema["step_num_col"], sensor_cols)

    dq_vals_check = [rs.get("mean_dq") for rs in run_signals.values()
                     if rs.get("mean_dq") is not None]
    dq_iqr_check  = _iqr(dq_vals_check) if len(dq_vals_check) >= 4 else 0.0
    if dq_iqr_check > 0:
        print(f"  DATA_QUALITY varies (IQR={dq_iqr_check:.2f}) -- PATH A highflier classification")
    else:
        print("  DATA_QUALITY is flat -- PATH B (fallback) highflier classification")

    # Flat abort set for tool-unaware step functions (run_id strings only).
    # A run aborted in any tool is excluded from step scoring entirely.
    aborted_flat = {rid for _, rid in aborted_all}

    step_info_si, run_order_si, run_tool_si, _ = score_steps(
        rows, schema["run_col"],
        schema["group_cols"][0] if schema.get("group_cols") else None,
        schema["step_num_col"], sensor_cols, aborted_flat,
        noisy_sensor_steps=schema.get("noisy_sensor_steps", set()))

    top_steps_si   = set(sorted(step_info_si, key=lambda s: -step_info_si[s]["mean_cv"])[:TOP_N_STEPS])

    # Build forced_step_sensors: for every anomalous run, determine which
    # (step, sensor) pairs are actually anomalous by rescoring per-step.
    # Only force a sensor into the trend chart for steps where its per-step
    # anomaly score is genuinely elevated for that specific anomalous run.
    forced_step_sensors = defaultdict(set)   # { step_number: {sensor, ...} }
    if schema.get("step_num_col"):
        for tool_id, rid in anom_set_all:
            psc      = per_sensor_sc_all.get(tool_id, {}).get(rid, {})
            baselines= baselines_full_all.get(tool_id, {})
            rrows_f  = tool_runs.get(tool_id, {}).get(rid, [])
            if not psc or not rrows_f:
                continue
            # Sensors with non-trivial composite score for this anomalous run
            flagged_sensors_run = {s for s, sc in psc.items() if sc > 0.0}
            # Collect all scored steps for this run
            steps_in_run = {str(r.get("_step_number", "")).strip()
                            for r in rrows_f
                            if str(r.get("_step_number", "")).strip()
                            and str(r.get("_step_number", "")).strip() not in SKIP_STEPS}
            for sn_f in steps_in_run:
                sf = {sn_f}
                for s in flagged_sensors_run:
                    bl = baselines.get(s)
                    if bl is None:
                        continue
                    # Score this sensor for this run restricted to this step
                    step_sc = compute_run_score(rrows_f, bl, s, step_filter=sf)
                    if step_sc > 0.0:
                        forced_step_sensors[sn_f].add(s)

    # Merge forced steps into top_steps_si
    top_steps_si.update(forced_step_sensors.keys())

    # Collect sensors for all top + forced steps
    all_si_sensors = sorted(
        {s for sn in top_steps_si for s, *_ in step_info_si.get(sn, {}).get("sensor_cvs", [])}
        | {s for sensors in forced_step_sensors.values() for s in sensors}
    )
    _tool_col_si = schema["group_cols"][0] if schema.get("group_cols") else None
    run_means_by_step = collect_run_means(
        rows, schema["run_col"], schema["step_num_col"],
        all_si_sensors, top_steps_si, aborted_flat,
        tool_col=_tool_col_si)

    # Per-run row-count lookup: ckey is "tool::rid"; look up (tool, rid) in row_counts_all.
    row_counts_flat = {}
    for ckey in run_order_si:
        tool_si = run_tool_si.get(ckey, "ALL")
        rid_si  = ckey.split("::", 1)[1] if "::" in ckey else ckey
        row_counts_flat[ckey] = row_counts_all.get((tool_si, rid_si), 0)

    # aborted_flat as composite keys for svg_trend_chart
    aborted_flat_ckeys = {f"{tool}::{rid}" for tool, rid in aborted_all}

    highfliers = detect_highfliers(
        run_means_by_step, run_order_si,
        run_signals, 0, median_step_rows)

    for sn, hf_step in highfliers.items():
        total = sum(len(f) for f in hf_step.values())
        if total:
            real = sum(1 for f in hf_step.values()
                       for h in f.values() if h["class"] == "real")
            print(f"  Step {sn}: {total} highflier(s) ({real} real, {total-real} artifact)")

    # ── Write step sections and close ────────────────────────────────────────
    print("\nWriting step importance sections (phase 2 complete) ...")
    _write_step_sections(fh, step_info_si, run_order_si, run_tool_si,
                          run_means_by_step, highfliers,
                          aborted_flat_ckeys, row_counts_flat,
                          forced_step_sensors=forced_step_sensors,
                          noisy_sensor_steps=schema.get("noisy_sensor_steps", set()))

    print(f"\nDone. Full report: {html_out}")


if __name__ == "__main__":
    main()
