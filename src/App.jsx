import { useState, useRef, useEffect, useCallback } from "react";
import { flushSync } from "react-dom";
import DatasetSlot from "./components/DatasetSlot";
import LogEntry from "./components/LogEntry";
import HypothesesPanel from "./components/HypothesesPanel";
import CoverageDock from "./components/CoverageDock";
import ReportDrawer from "./components/ReportDrawer";
import { setGroupMappings, uploadDegDataset, deleteDegDataset } from "./api";
import { THEMES, FONT_SANS, RADII, SHADOW, cssVars, ACCENTS, applyAccent } from "./theme";

function makeStyles(t) {
  return `
  @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;500;600&display=swap');
  *{box-sizing:border-box;margin:0;padding:0}
  html,body,#root{background:${t.appBg};width:100%;height:100%;overflow:hidden}
  ::-webkit-scrollbar{width:5px}
  ::-webkit-scrollbar-track{background:transparent}
  ::-webkit-scrollbar-thumb{background:${t.border};border-radius:3px}
  ::-webkit-scrollbar-thumb:hover{background:${t.accent}55}
  @keyframes pulse{0%,100%{opacity:1}50%{opacity:.35}}
  @keyframes si{from{opacity:0;transform:translateY(3px)}to{opacity:1;transform:translateY(0)}}
  @keyframes dots{0%,100%{content:''}33%{content:'.'}66%{content:'..'}99%{content:'...'}}
  .thinking-indicator::after{content:'';animation:dots 1.2s steps(1) infinite}
  @keyframes spin{from{transform:rotate(0deg)}to{transform:rotate(360deg)}}
  .spinner{width:14px;height:14px;border:2px solid ${t.border};border-top-color:${t.accent};border-radius:50%;animation:spin 0.75s linear infinite;flex-shrink:0}
  @keyframes thb{0%,100%{transform:scaleY(.35)}50%{transform:scaleY(1)}}
  .thrb{display:inline-flex;align-items:flex-end;gap:2px;height:11px}
  .thrb i{width:2px;height:100%;background:currentColor;border-radius:1px;transform-origin:bottom;animation:thb 1s ease-in-out infinite}
  .ent{animation:si .18s ease}
  .blink{animation:pulse 1.4s infinite}
  .lockzone{border:0;margin:0;padding:0;min-width:0}
  .lockzone:disabled{opacity:.55}
  .btn{
    background:${t.cardBg};border:1px solid ${t.border};color:${t.textPrimary};
    font-family:inherit;font-size:13px;padding:8px 14px;cursor:pointer;
    transition:background .12s,border-color .12s,color .12s;
    width:100%;border-radius:${RADII.md}px;font-weight:500;
  }
  .btn:hover{background:${t.surface2};border-color:${t.borderStrong}}
  .btn:disabled{opacity:.4;cursor:not-allowed}
  .bsm{padding:5px 10px;width:auto;font-size:12px}
  .bdng{border-color:${t.warning}55;color:${t.warning};background:transparent}
  .bdng:hover{background:${t.dangerHoverBg};border-color:${t.warning}}
  .slot{
    border:1px solid ${t.border};padding:11px;margin-bottom:8px;
    background:${t.cardBg};border-radius:${RADII.md}px;transition:border-color .12s,box-shadow .12s;
  }
  .slot:hover{border-color:${t.borderStrong};box-shadow:var(--shadow-sm)}
  .slot.ok{border-color:${t.confirmedBd}}
  .uz{
    border:1px dashed ${t.borderStrong};padding:9px 12px;text-align:center;cursor:pointer;
    transition:all .12s;background:transparent;display:flex;align-items:center;
    justify-content:center;gap:6px;margin-bottom:6px;font-size:12px;
    color:${t.textMuted};border-radius:${RADII.md}px;
  }
  .uz:hover{border-color:${t.accent};background:${t.accentSoft};color:${t.accent}}
  .uz.ok{border-color:${t.confirmedBd};color:${t.confirmed}}
  .tag{display:inline-block;padding:2px 7px;border-radius:${RADII.sm}px;font-size:11px;font-weight:600}
  input[type=text],input[type=number],select{
    background:${t.cardBg};border:1px solid ${t.border};color:${t.textPrimary};
    padding:8px 10px;font-size:12.5px;font-family:inherit;width:100%;
    border-radius:${RADII.md}px;transition:border-color .12s;
  }
  input[type=text]:focus,input[type=number]:focus,select:focus{outline:none;border-color:${t.borderStrong}}
  .sec{
    font-size:10.5px;color:${t.textMuted};letter-spacing:0.06em;
    margin:18px 0 9px;font-weight:600;text-transform:uppercase;
    display:flex;align-items:center;gap:8px;
  }
  .sec:first-child{margin-top:6px}
  `;
}

function SettingRow({ label, children }) {
  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 7 }}>
      <span style={{ fontSize: 10, fontWeight: 700, letterSpacing: "0.06em", textTransform: "uppercase", color: "var(--text-3)" }}>{label}</span>
      {children}
    </div>
  );
}

function Seg({ options, value, onChange, t }) {
  return (
    <div style={{ display: "flex", gap: 2, padding: 3, background: t.appBg, border: `1px solid ${t.border}`, borderRadius: RADII.md }}>
      {options.map(({ k, l }) => (
        <button key={k} onClick={() => onChange(k)}
          style={{ flex: 1, padding: "5px 6px", fontSize: 11.5, fontFamily: "inherit", cursor: "pointer", borderRadius: RADII.sm,
            background: value === k ? t.cardBg : "transparent",
            border: value === k ? `1px solid ${t.accent}40` : "1px solid transparent",
            color: value === k ? t.accent : t.textMuted, fontWeight: value === k ? 600 : 400 }}>
          {l}
        </button>
      ))}
    </div>
  );
}

export default function App() {
  const [colorMode, setColorMode] = useState(() => localStorage.getItem("ta_theme") || "dark");
  const [accent,    setAccent]    = useState(() => localStorage.getItem("ta_accent") || "graphite");
  const [density,   setDensity]   = useState(() => localStorage.getItem("ta_density") || "comfortable");
  const [settingsOpen, setSettingsOpen] = useState(false);
  const [reportOpen, setReportOpen] = useState(false);
  const t = applyAccent(THEMES[colorMode], colorMode, accent);

  useEffect(() => { localStorage.setItem("ta_theme", colorMode); }, [colorMode]);
  useEffect(() => { localStorage.setItem("ta_accent", accent); }, [accent]);
  useEffect(() => { localStorage.setItem("ta_density", density); }, [density]);

  const [slots, setSlots] = useState([
    { id: 0, exprFile: null, metaFile: null, name: "Dataset 1" },
  ]);
  const [loaded,        setLoaded]        = useState([]);
  const [groupMap,      setGroupMap]      = useState({});
  const [phase,         setPhase]         = useState("upload");
  const [log,           setLog]           = useState([]);
  const [hypotheses,    setHypotheses]    = useState([]);
  const [step,          setStep]          = useState(0);
  const [maxHypotheses, setMaxHypotheses] = useState(3);
  const [agentMode,     setAgentMode]     = useState("reproduce");
  const [piModel,       setPiModel]       = useState("claude-opus-4-8");
  const [currentStatus, setCurrentStatus] = useState("");
  const [streamingText, setStreamingText] = useState("");
  const [mappingGroups, setMappingGroups] = useState([]);
  const [mappingsOpen,  setMappingsOpen]  = useState(false);
  const [degDatasets,   setDegDatasets]   = useState([]);
  const [degFile,       setDegFile]       = useState(null);
  const [degGroupA,     setDegGroupA]     = useState("");
  const [degGroupB,     setDegGroupB]     = useState("");
  const [degUploading,  setDegUploading]  = useState(false);
  const [degStatus,     setDegStatus]     = useState("");
  const [runCost,       setRunCost]       = useState(null);
  const [priorCount,    setPriorCount]    = useState(0);
  const logEnd    = useRef(null);
  const scrollRef = useRef(null);
  const abortRef  = useRef(null);
  const [atBottom, setAtBottom] = useState(true);

  // Auto-scroll to the newest entry only while the user is pinned to the bottom,
  // so reading back through the log isn't yanked away as new steps stream in.
  useEffect(() => {
    const el = scrollRef.current;
    if (el && atBottom) el.scrollTop = el.scrollHeight;
  }, [log, streamingText, currentStatus, atBottom]);

  const onLogScroll = () => {
    const el = scrollRef.current;
    if (!el) return;
    setAtBottom(el.scrollHeight - el.scrollTop - el.clientHeight < 80);
  };
  const jumpToLatest = () => {
    const el = scrollRef.current;
    if (el) { el.scrollTop = el.scrollHeight; setAtBottom(true); }
  };
  const addLog = useCallback(e => setLog(prev => [...prev, { ...e, id: Date.now() + Math.random() }]), []);

  const addSlot    = () => setSlots(p => [...p, { id: Date.now(), exprFile: null, metaFile: null, name: `Dataset ${p.length + 1}` }]);
  const removeSlot = id => setSlots(p => p.filter(s => s.id !== id));
  const updSlot    = (id, k, v) => setSlots(p => p.map(s => s.id === id ? { ...s, [k]: v } : s));

  const postMappings = useCallback((groups) => {
    const mappings = {};
    for (const { canonical, aliases } of groups) {
      if (canonical.trim() && aliases.size > 0) mappings[canonical.trim()] = [...aliases];
    }
    setGroupMappings(mappings).catch(() => {});
  }, []);

  const addMappingGroup = () => {
    const next = [...mappingGroups, { canonical: "", aliases: new Set() }];
    setMappingGroups(next);
    postMappings(next);
  };

  const updateMappingCanonical = (idx, value) => {
    const next = mappingGroups.map((mg, i) => i === idx ? { ...mg, canonical: value } : mg);
    setMappingGroups(next);
    postMappings(next);
  };

  const toggleAlias = (idx, alias, checked) => {
    const next = mappingGroups.map((mg, i) => {
      if (i !== idx) return mg;
      const aliases = new Set(mg.aliases);
      if (checked) aliases.add(alias); else aliases.delete(alias);
      return { ...mg, aliases };
    });
    setMappingGroups(next);
    postMappings(next);
  };

  const removeMappingGroup = (idx) => {
    const next = mappingGroups.filter((_, i) => i !== idx);
    setMappingGroups(next);
    postMappings(next);
  };

  const loadAll = async () => {
    const newLoaded = [];
    for (const slot of slots) {
      if (!slot.exprFile || !slot.metaFile) continue;
      const fd = new FormData();
      fd.append("expr_file", slot.exprFile);
      fd.append("meta_file", slot.metaFile);
      fd.append("name", slot.name);
      const existingId = loaded.find(d => d.name === slot.name)?.id;
      if (existingId) fd.append("dataset_id", existingId);
      try {
        const res = await fetch("/api/datasets", { method: "POST", body: fd });
        if (!res.ok) { addLog({ type: "error", text: `Upload ${slot.name}: ${res.statusText}` }); continue; }
        newLoaded.push(await res.json());
      } catch (e) {
        addLog({ type: "error", text: `Upload ${slot.name}: ${e.message}` });
      }
    }
    setLoaded(newLoaded);
    const m = {};
    newLoaded.forEach(d => { m[d.id] = d.group_col; });
    setGroupMap(m);
    setMappingsOpen(newLoaded.length >= 2);
  };

  const removeDeg = (name) => {
    setDegDatasets(prev => prev.filter(x => x.name !== name));
    deleteDegDataset(name).catch(() => {});
  };

  const uploadDeg = async () => {
    if (!degFile || !degGroupA.trim() || !degGroupB.trim()) return;
    // Pick the lowest unused index so removing then re-adding a table can't
    // collide with an existing name and silently overwrite it server-side.
    const used = new Set(
      degDatasets.map(d => parseInt((d.name.match(/\d+/) || [])[0], 10)).filter(Number.isFinite)
    );
    let n = 1;
    while (used.has(n)) n++;
    const degName = `DEG ${n}`;
    setDegUploading(true);
    setDegStatus("");
    try {
      const res = await uploadDegDataset(degFile, degName, degGroupA.trim(), degGroupB.trim());
      const errMsg = res.error || res.detail;
      if (errMsg) { setDegStatus(`Error: ${errMsg}`); return; }
      setDegDatasets(prev => [...prev.filter(d => d.name !== res.name), res]);
      setDegStatus(`Uploaded: ${res.n_genes} genes (${res.groupA} vs ${res.groupB})`);
      setMappingsOpen(true);
      setDegFile(null);
      setDegGroupA("");
      setDegGroupB("");
    } catch (e) {
      setDegStatus(`Error: ${e.message}`);
    } finally {
      setDegUploading(false);
    }
  };

  const runAgent = async () => {
    if (!loaded.length && !degDatasets.length) return;
    setPhase("running"); setLog([]); setStep(0); setHypotheses([]); setCurrentStatus("Running pre-analysis..."); setStreamingText(""); setRunCost(0); setPriorCount(0);

    const controller = new AbortController();
    abortRef.current = controller;
    let currentStep = 0;

    try {
      const res = await fetch("/api/run", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ dataset_ids: loaded.map(d => d.id), group_cols: groupMap, max_hypotheses: maxHypotheses, mode: agentMode, model: piModel }),
        signal: controller.signal,
      });

      if (!res.ok) {
        const err = await res.json().catch(() => ({ detail: res.statusText }));
        addLog({ type: "error", text: `Start: ${err.detail || res.statusText}` });
        setPhase("done");
        return;
      }

      const reader = res.body.getReader();
      const decoder = new TextDecoder();
      let buf = "";

      while (true) {
        const { done, value } = await reader.read();
        if (done) break;
        buf += decoder.decode(value, { stream: true });
        const lines = buf.split("\n");
        buf = lines.pop();
        for (const line of lines) {
          if (!line.startsWith("data: ")) continue;
          try {
            const entry = JSON.parse(line.slice(6));
            if (entry.type === "stream_end") { reader.cancel(); break; }
            if (entry.type === "thinking")       { setStep(++currentStep); setCurrentStatus(entry.text); setStreamingText(""); addLog(entry); continue; }
            if (entry.type === "thought_stream") { flushSync(() => setStreamingText(prev => prev + entry.delta)); continue; }
            if (entry.type === "thought")        { flushSync(() => setStreamingText("")); }
            if (entry.type === "usage")          { setRunCost(entry.total_cost_usd); continue; }
            if (entry.type === "prior_knowledge") { setPriorCount(entry.count || 0); continue; }
            if (entry.type === "hypothesis_propose") setHypotheses(prev => [...prev, entry.hypothesis]);
            if (entry.type === "hypothesis_eval")    setHypotheses(prev => prev.map(h => h.id === entry.hypothesis.id ? entry.hypothesis : h));
            addLog(entry);
            await new Promise(r => setTimeout(r, 80));
          } catch { /* ignore malformed SSE lines */ }
        }
      }
    } catch (e) {
      if (e.name !== "AbortError") addLog({ type: "error", text: `Stream: ${e.message}` });
    }
    setCurrentStatus("");
    setStreamingText("");
    setPhase("done");
  };

  const hasData = loaded.length > 0 || degDatasets.length > 0;
  const locked = phase === "running";   // setup is read-only while a run is in flight

  return (
    <div style={{ ...cssVars(colorMode), "--accent": t.accent, "--accent-hover": t.accentHover, "--accent-soft": t.accentSoft, "--accent-text-on": t.accentTextOn, minHeight: "100vh", background: t.appBg, fontFamily: FONT_SANS, fontSize: density === "compact" ? 12.5 : 13.5, color: t.textPrimary }}>
      <style>{makeStyles(t)}</style>

      {/* Header */}
      <div style={{ height: 52, padding: "0 16px", display: "flex", alignItems: "center", gap: 14, background: t.sidebarBg, flexShrink: 0, borderBottom: `1px solid ${t.border}`, zIndex: 30 }}>

        {/* brand */}
        <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
          <div style={{ width: 31, height: 31, borderRadius: 8, background: t.accent, color: t.accentTextOn, display: "grid", placeItems: "center", flexShrink: 0, boxShadow: SHADOW[colorMode].sm }}>
            <svg viewBox="0 0 24 24" fill="none" style={{ width: 20, height: 20 }}>
              <circle cx="12" cy="8.7" r="6.05" fill="currentColor" fillOpacity="0.5" />
              <circle cx="8.4" cy="15" r="6.05" fill="currentColor" fillOpacity="0.5" />
              <circle cx="15.6" cy="15" r="6.05" fill="currentColor" fillOpacity="0.5" />
            </svg>
          </div>
          <div style={{ display: "flex", flexDirection: "column", lineHeight: 1.04 }}>
            <span style={{ fontWeight: 600, letterSpacing: "-0.02em", fontSize: 15.5, color: t.textPrimary }}>Transcriptomic Agent</span>
            <span style={{ fontFamily: "'IBM Plex Mono',ui-monospace,monospace", fontSize: 8.5, letterSpacing: "0.16em", textTransform: "uppercase", color: t.textMuted, marginTop: 1 }}>Discovery engine</span>
          </div>
        </div>

        <div style={{ width: 1, height: 22, background: t.border }} />

        {/* run status pill */}
        {(phase === "running" || phase === "done") && (
          <div style={{ display: "flex", alignItems: "center", gap: 8, whiteSpace: "nowrap", fontSize: 12, color: t.textSecondary, padding: "4px 11px 4px 9px", border: `1px solid ${t.border}`, borderRadius: 99, background: t.surface2 }}>
            <span className={phase === "running" ? "blink" : ""} style={{ width: 7, height: 7, borderRadius: 99, background: phase === "running" ? t.accent : t.confirmed, flexShrink: 0 }} />
            {phase === "running"
              ? <>Running · {hypotheses.filter(h => h.seeded_by === "llm" && h.status !== "pending").length}/{maxHypotheses} off-grid</>
              : <>Completed · {hypotheses.filter(h => ["confirmed", "uncertain", "rejected"].includes(h.status)).length} adjudicated</>}
          </div>
        )}

        <div style={{ flex: 1 }} />

        {priorCount > 0 && (
          <div title="Confirmed findings recalled from prior runs on this data" style={{ display: "flex", alignItems: "center", gap: 6, padding: "4px 11px", fontSize: 12, color: t.textSecondary, background: t.surface2, border: `1px solid ${t.border}`, borderRadius: 99 }}>
            <span style={{ fontSize: 12 }}>🧠</span>
            Remembers <b style={{ color: t.textPrimary, fontWeight: 600 }}>{priorCount}</b> prior finding{priorCount === 1 ? "" : "s"}
          </div>
        )}
        {runCost !== null && (
          <div title={`Estimated API cost (${piModel})`} style={{ display: "flex", alignItems: "center", gap: 5, padding: "4px 11px", background: t.surface2, border: `1px solid ${t.border}`, borderRadius: 99 }}>
            <span style={{ fontSize: 11, color: t.textMuted, fontFamily: "'IBM Plex Mono',ui-monospace,monospace" }}>$</span>
            <span style={{ fontSize: 12, color: t.textSecondary, fontFamily: "'IBM Plex Mono',ui-monospace,monospace", fontWeight: 500 }}>
              {runCost < 0.01 ? runCost.toFixed(4) : runCost.toFixed(3)}
            </span>
          </div>
        )}
        {hypotheses.length > 0 && (
          <button
            title="Open findings report"
            onClick={() => setReportOpen(true)}
            style={{ height: 32, padding: "0 11px", display: "inline-flex", alignItems: "center", gap: 6, background: "transparent", border: `1px solid ${t.accent}40`, color: t.accent, fontSize: 12, borderRadius: RADII.md, cursor: "pointer", fontFamily: "inherit", transition: "all .15s" }}
            onMouseEnter={e => { e.currentTarget.style.background = t.startHoverBg; }}
            onMouseLeave={e => { e.currentTarget.style.background = "transparent"; }}
          >
            📄 Report
          </button>
        )}
        <div style={{ position: "relative" }}>
          <button
            title="Settings"
            onClick={() => setSettingsOpen(o => !o)}
            style={{ height: 32, padding: "0 11px", display: "inline-flex", alignItems: "center", gap: 6, background: settingsOpen ? t.surface2 : "transparent", border: `1px solid ${t.border}`, color: settingsOpen ? t.textPrimary : t.textSecondary, fontSize: 12, borderRadius: RADII.md, cursor: "pointer", fontFamily: "inherit", transition: "all .15s" }}
            onMouseEnter={e => { e.currentTarget.style.background = t.surface2; e.currentTarget.style.color = t.textPrimary; }}
            onMouseLeave={e => { if (!settingsOpen) { e.currentTarget.style.background = "transparent"; e.currentTarget.style.color = t.textSecondary; } }}
          >
            ⚙ Settings
          </button>

          {settingsOpen && (
            <>
              <div onClick={() => setSettingsOpen(false)} style={{ position: "fixed", inset: 0, zIndex: 40 }} />
              <div style={{ position: "absolute", top: 40, right: 0, zIndex: 41, width: 232, padding: 14, background: t.cardBg, border: `1px solid ${t.border}`, borderRadius: RADII.lg, boxShadow: SHADOW[colorMode].lg, display: "flex", flexDirection: "column", gap: 14 }}>
                <SettingRow label="Theme">
                  <Seg options={[{ k: "light", l: "Light" }, { k: "dark", l: "Dark" }]} value={colorMode} onChange={setColorMode} t={t} />
                </SettingRow>
                <SettingRow label="Accent">
                  <div style={{ display: "flex", gap: 8 }}>
                    {Object.entries(ACCENTS).map(([key, a]) => (
                      <button key={key} title={a.label} onClick={() => setAccent(key)}
                        style={{ width: 30, height: 30, borderRadius: 8, background: a.swatch, cursor: "pointer", padding: 0, border: accent === key ? `2px solid ${t.textPrimary}` : `2px solid transparent`, boxShadow: accent === key ? `0 0 0 2px ${t.cardBg} inset` : "none" }} />
                    ))}
                  </div>
                </SettingRow>
                <SettingRow label="Density">
                  <Seg options={[{ k: "compact", l: "Compact" }, { k: "comfortable", l: "Comfortable" }]} value={density} onChange={setDensity} t={t} />
                </SettingRow>
              </div>
            </>
          )}
        </div>
      </div>

      <div style={{ display: "flex", height: "calc(100vh - 52px)" }}>

        {/* LEFT PANEL */}
        <div style={{ width: 288, borderRight: `1px solid ${t.border}`, padding: "8px 14px 14px", overflowY: "auto", flexShrink: 0, background: t.sidebarBg }}>

          {locked && (
            <div style={{ display: "flex", alignItems: "center", gap: 7, margin: "4px 0 10px", padding: "6px 10px", fontSize: 11.5, color: t.textSecondary, background: t.surface2, border: `1px solid ${t.border}`, borderRadius: RADII.sm }}>
              <span style={{ fontSize: 12 }}>🔒</span> Setup locked while the run is in flight
            </div>
          )}

          <fieldset className="lockzone" disabled={locked}>
          <div className="sec">Datasets</div>

          {slots.map(slot => (
            <DatasetSlot key={slot.id} slot={slot} canRemove={slots.length > 1} theme={t}
              onUpdate={(k, v) => updSlot(slot.id, k, v)} onRemove={() => removeSlot(slot.id)} />
          ))}

          <button className="btn bsm" style={{ marginBottom: 6, width: "100%" }} onClick={addSlot}>+ Add dataset</button>
          <button className="btn" style={{ marginBottom: 4, borderColor: `${t.accent}30`, color: t.accent }} onClick={loadAll} disabled={!slots.some(s => s.exprFile && s.metaFile)}>
            Load datasets
          </button>

          {/* DEG TABLE UPLOAD */}
          <div className="sec">DEG Datasets</div>
          <div style={{ marginBottom: 4 }}>
            <label className="uz" style={{ marginBottom: 6 }}>
              {degFile ? degFile.name : "Upload DEG table (.csv)"}
              <input type="file" accept=".csv" style={{ display: "none" }} onChange={e => setDegFile(e.target.files[0] || null)} />
            </label>
            <input type="text" value={degGroupA} onChange={e => setDegGroupA(e.target.value)} placeholder="Group A" style={{ marginBottom: 5 }} />
            <input type="text" value={degGroupB} onChange={e => setDegGroupB(e.target.value)} placeholder="Group B" style={{ marginBottom: 6 }} />
            <button className="btn bsm" style={{ width: "100%", marginBottom: 6 }}
              onClick={uploadDeg} disabled={!degFile || !degGroupA.trim() || !degGroupB.trim() || degUploading}>
              {degUploading ? "Uploading..." : "Upload DEG table"}
            </button>
            {degStatus && (
              <div style={{ fontSize: 12, marginBottom: 8, padding: "5px 9px", borderRadius: RADII.sm,
                color: degStatus.startsWith("Error") ? t.warning : t.confirmed,
                background: degStatus.startsWith("Error") ? t.warningSoft : t.confirmedSoft,
                border: `1px solid ${degStatus.startsWith("Error") ? `${t.warning}40` : t.confirmedBd}` }}>
                {degStatus}
              </div>
            )}
            {degDatasets.map(d => (
              <div key={d.name} className="slot ok" style={{ marginBottom: 8 }}>
                <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 6 }}>
                  <span style={{ fontSize: 13, color: t.accent, fontWeight: 600 }}>{d.name}</span>
                  <button className="btn bsm bdng" style={{ padding: "2px 8px", fontSize: 11 }}
                    onClick={() => removeDeg(d.name)}>✕</button>
                </div>
                {(d.comparisons || []).map((c, i) => (
                  <div key={i} style={{ fontSize: 12, color: t.textSecondary, lineHeight: 1.8, fontFamily: "'IBM Plex Mono',ui-monospace,monospace" }}>
                    {c.groupA} <span style={{ color: t.textMuted }}>vs</span> {c.groupB}
                    <span style={{ color: t.textMuted, marginLeft: 6 }}>{c.n_genes} genes</span>
                  </div>
                ))}
              </div>
            ))}
          </div>

          </fieldset>

          {hasData && <>
            <fieldset className="lockzone" disabled={locked}>
            <div className="sec">Group Columns</div>
            {loaded.map(ds => (
              <div key={ds.id} style={{ marginBottom: 14 }}>
                <div style={{ fontSize: 12, color: t.textSecondary, marginBottom: 5, fontWeight: 600 }}>{ds.name}</div>
                <select value={groupMap[ds.id] || ds.group_col} onChange={async e => {
                    const newCol = e.target.value;
                    setGroupMap(prev => ({ ...prev, [ds.id]: newCol }));
                    try {
                      const r = await fetch(`/api/datasets/${ds.id}/group_col`, {
                        method: "PATCH",
                        headers: { "Content-Type": "application/json" },
                        body: JSON.stringify({ group_col: newCol }),
                      });
                      if (r.ok) {
                        const data = await r.json();
                        setLoaded(prev => prev.map(d => d.id === ds.id ? { ...d, group_col: newCol, groups: data.groups, group_counts: data.group_counts } : d));
                      }
                    } catch {}
                  }}>
                  {(ds.group_col_candidates || (ds.group_cols || []).map(c => ({ col: c, unique_values: [] }))).map(cand => (
                    <option key={cand.col} value={cand.col}>
                      {cand.col}: {cand.unique_values.slice(0, 3).join(", ")}
                    </option>
                  ))}
                </select>
                <div style={{ fontSize: 11, color: t.textMuted, marginTop: 5, lineHeight: 1.9, fontFamily: "'IBM Plex Mono',ui-monospace,monospace" }}>
                  {ds.groups.map(g => (
                    <div key={g} style={{ paddingLeft: 2, display: "flex", justifyContent: "space-between", gap: 8 }}>
                      <span>{g}</span>
                      {ds.group_counts?.[g] != null && <span style={{ color: t.textFaint }}>n={ds.group_counts[g]}</span>}
                    </div>
                  ))}
                </div>
                <div style={{ fontSize: 11, color: t.textMuted, marginTop: 4, opacity: 0.7 }}>
                  {ds.gene_count} genes · {ds.sample_count} samples
                </div>
              </div>
            ))}

            {/* GROUP MAPPINGS */}
            {(() => {
              const degGroups = degDatasets.flatMap(d => {
                const fromComparisons = (d.comparisons || []).flatMap(c => [c.groupA, c.groupB]);
                const topLevel = [d.groupA, d.groupB].filter(Boolean);
                return [...fromComparisons, ...topLevel];
              });
              const allGroups = [...new Set([...loaded.flatMap(ds => ds.groups), ...degGroups])];
              return (
                <>
                  <div className="sec" style={{ cursor: "pointer", userSelect: "none" }}
                    onClick={() => setMappingsOpen(p => !p)}>
                    <span>Group Mappings</span>
                    <span style={{ color: t.accent, fontSize: 11, marginLeft: -4 }}>{mappingsOpen ? "▾" : "▸"}</span>
                  </div>
                  {mappingsOpen && (
                    <div style={{ marginBottom: 10 }}>
                      {mappingGroups.map((mg, idx) => (
                        <div key={idx} style={{ marginBottom: 8, padding: "9px 10px", border: `1px solid ${t.border}`, background: t.cardBg, borderRadius: 7 }}>
                          <div style={{ display: "flex", gap: 5, marginBottom: 7 }}>
                            <input type="text" value={mg.canonical} placeholder="Canonical name"
                              onChange={e => updateMappingCanonical(idx, e.target.value)} style={{ flex: 1 }} />
                            <button className="btn bsm bdng" onClick={() => removeMappingGroup(idx)}>✕</button>
                          </div>
                          <div style={{ display: "flex", flexWrap: "wrap", gap: 5, marginBottom: 5 }}>
                            {allGroups.map(g => (
                              <label key={g} style={{ fontSize: 12, color: mg.aliases.has(g) ? t.accent : t.textMuted, cursor: "pointer", display: "flex", alignItems: "center", gap: 4 }}>
                                <input type="checkbox" checked={mg.aliases.has(g)} onChange={e => toggleAlias(idx, g, e.target.checked)} />
                                {g}
                              </label>
                            ))}
                          </div>
                          {mg.canonical && (
                            <div style={{ fontSize: 11, color: t.textMuted, lineHeight: 1.6, fontFamily: "'IBM Plex Mono',ui-monospace,monospace" }}>
                              "{mg.canonical}" ← {allGroups.map(g => (
                                <span key={g} style={{ marginRight: 6, color: mg.aliases.has(g) ? t.accent : t.textMuted }}>
                                  {g} {mg.aliases.has(g) ? "✓" : "✗"}
                                </span>
                              ))}
                            </div>
                          )}
                        </div>
                      ))}
                      <button className="btn bsm" style={{ width: "100%", marginBottom: 6 }} onClick={addMappingGroup}>
                        + Add mapping group
                      </button>
                    </div>
                  )}
                </>
              );
            })()}

            <div className="sec">Mode</div>
            <div style={{ display: "flex", background: t.appBg, borderRadius: 7, border: `1px solid ${t.border}`, padding: 3, marginBottom: 12, gap: 2 }}>
              {[
                { key: "reproduce", label: "Reproduce", sub: "deterministic" },
                { key: "explore",   label: "Explore",   sub: "creative"       },
              ].map(({ key, label, sub }) => (
                <button key={key} onClick={() => setAgentMode(key)}
                  style={{
                    flex: 1,
                    background: agentMode === key ? t.cardBg : "transparent",
                    border: agentMode === key ? `1px solid ${t.accent}30` : "1px solid transparent",
                    color: agentMode === key ? t.accent : t.textMuted,
                    padding: "6px 8px", cursor: "pointer", borderRadius: 5, transition: "all .15s",
                    fontFamily: "inherit", fontSize: 12, fontWeight: agentMode === key ? 600 : 400,
                  }}>
                  {label}
                  <div style={{ fontSize: 10, color: agentMode === key ? t.accent : t.textMuted, marginTop: 1, opacity: agentMode === key ? 0.7 : 0.5 }}>{sub}</div>
                </button>
              ))}
            </div>

            <div className="sec">PI Model</div>
            <select value={piModel} onChange={e => setPiModel(e.target.value)} style={{ marginBottom: 12 }}>
              <option value="claude-opus-4-8">Claude Opus 4.8 — most capable</option>
              <option value="claude-sonnet-4-6">Claude Sonnet 4.6 — balanced</option>
              <option value="claude-haiku-4-5">Claude Haiku 4.5 — fast & cheap</option>
            </select>

            <div className="sec">Off-grid hypotheses</div>
            <input type="number" value={maxHypotheses} min={0} max={30}
              onChange={e => setMaxHypotheses(Math.max(0, parseInt(e.target.value) || 0))}
              style={{ marginBottom: 6 }} />
            <div style={{ fontSize: 11, color: t.textMuted, marginBottom: 14, lineHeight: 1.5 }}>
              Budget for agent-proposed (off-grid) hypotheses. The deterministic Layer 1 — pairwise
              comparisons and the coverage grid — always runs in full regardless and is not counted here.
              Set to 0 for a Layer-1-only run.
            </div>
            </fieldset>

            <button
              style={{
                width: "100%", padding: "10px 14px", border: "1px solid transparent", borderRadius: 7,
                background: phase === "running" ? t.cardBg : t.elevatedBg,
                borderColor: phase === "running" ? t.border : `${t.accent}40`,
                color: phase === "running" ? t.textSecondary : t.accent,
                fontFamily: "inherit", fontSize: 13, fontWeight: 600,
                cursor: "pointer", transition: "all .15s", letterSpacing: 0.2,
              }}
              onMouseEnter={e => { if (phase !== "running") { e.currentTarget.style.background = t.startHoverBg; e.currentTarget.style.borderColor = `${t.accentHover}60`; e.currentTarget.style.color = t.accentHover; } }}
              onMouseLeave={e => { if (phase !== "running") { e.currentTarget.style.background = t.elevatedBg; e.currentTarget.style.borderColor = `${t.accent}40`; e.currentTarget.style.color = t.accent; } }}
              onClick={phase === "running" ? () => abortRef.current?.abort() : runAgent}>
              {phase === "running" ? "Stop" : "Start Agent"}
            </button>
          </>}
        </div>

        {/* LOG PANEL */}
        <div style={{ flex: 1, display: "flex", flexDirection: "column", overflow: "hidden", position: "relative" }}>

          {currentStatus && (
            <div style={{ flexShrink: 0, display: "flex", alignItems: "center", gap: 10, padding: "9px 24px", borderBottom: `1px solid ${t.border}`, background: t.sidebarBg }}>
              <div className="spinner" />
              <span className="thinking-indicator" style={{ fontSize: 13, color: t.accent }}>{currentStatus}</span>
              <div style={{ marginLeft: "auto", display: "flex", alignItems: "center", gap: 16 }}>
                {phase === "running" && (
                  <span style={{ display: "inline-flex", alignItems: "center", gap: 6, fontSize: 11, color: t.textMuted }}>
                    <span className="thrb" style={{ color: t.accent }}>
                      <i style={{ animationDelay: "0s" }} /><i style={{ animationDelay: ".15s" }} /><i style={{ animationDelay: ".3s" }} /><i style={{ animationDelay: ".45s" }} />
                    </span>
                    streaming
                  </span>
                )}
                {hypotheses.length > 0 && <span style={{ fontSize: 12, color: t.textMuted }}>{hypotheses.filter(h => h.seeded_by === "llm" && h.status !== "pending").length}/{maxHypotheses} off-grid evaluated</span>}
              </div>
            </div>
          )}

          <div ref={scrollRef} onScroll={onLogScroll} style={{ flex: 1, overflowY: "auto", padding: "24px 28px" }}>
            {log.length === 0 && !currentStatus && (
              <div style={{ textAlign: "center", marginTop: "26vh" }}>
                <div style={{ width: 52, height: 52, margin: "0 auto 20px", background: `${t.accent}10`, border: `1px solid ${t.accent}20`, borderRadius: 14, display: "flex", alignItems: "center", justifyContent: "center", fontSize: 24, color: `${t.accent}50` }}>◈</div>
                <div style={{ fontSize: 15, color: t.textSecondary, fontWeight: 600, marginBottom: 8 }}>Ready to analyze</div>
                <div style={{ fontSize: 13, color: t.textMuted }}>Load datasets from the left panel, then start the agent!</div>
              </div>
            )}
            {log.length === 0 && currentStatus && (
              <div style={{ textAlign: "center", marginTop: "30vh" }}>
                <div style={{ display: "flex", justifyContent: "center", marginBottom: 20 }}>
                  <div style={{ width: 32, height: 32, border: `2px solid ${t.border}`, borderTopColor: t.accent, borderRadius: "50%", animation: "spin 0.75s linear infinite" }} />
                </div>
                <div style={{ fontSize: 14, color: t.accent, marginBottom: 6 }}>{currentStatus}</div>
                <div style={{ fontSize: 12, color: t.textMuted }}>This may take a few seconds...</div>
              </div>
            )}
            {log.map(e => <LogEntry key={e.id} entry={e} theme={t} />)}
            {streamingText && (
              <div className="ent" style={{ marginBottom: 12, paddingLeft: 16, borderLeft: `2px solid ${t.border}` }}>
                <div style={{ fontSize: 14, color: t.textPrimary, lineHeight: 1.75 }}>
                  {streamingText}<span className="blink" style={{ color: t.accent }}>▋</span>
                </div>
              </div>
            )}
            <div ref={logEnd} />
          </div>

          {!atBottom && log.length > 0 && (
            <button
              onClick={jumpToLatest}
              style={{ position: "absolute", left: "50%", transform: "translateX(-50%)", bottom: 16, display: "inline-flex", alignItems: "center", gap: 6, padding: "6px 13px", fontSize: 12, fontWeight: 500, color: t.accentTextOn, background: t.accent, border: "none", borderRadius: 99, cursor: "pointer", fontFamily: "inherit", boxShadow: SHADOW[colorMode].md, zIndex: 5 }}
            >
              ↓ Jump to latest
            </button>
          )}

          {(phase === "running" || hypotheses.length > 0) && <CoverageDock hypotheses={hypotheses} theme={t} />}
        </div>

        {/* HYPOTHESIS PANEL */}
        {(phase === "running" || hypotheses.length > 0) && (
          <div style={{ width: 300, borderLeft: `1px solid ${t.border}`, flexShrink: 0, background: t.sidebarBg }}>
            <HypothesesPanel hypotheses={hypotheses} maxHypotheses={maxHypotheses} theme={t} />
          </div>
        )}
      </div>

      {reportOpen && (
        <ReportDrawer hypotheses={hypotheses} datasets={loaded} colorMode={colorMode} theme={t} onClose={() => setReportOpen(false)} />
      )}
    </div>
  );
}
