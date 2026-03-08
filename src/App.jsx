import { useState, useRef, useEffect, useCallback } from "react";
import { flushSync } from "react-dom";
import DatasetSlot from "./components/DatasetSlot";
import LogEntry from "./components/LogEntry";

const STYLES = `
  @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;600&family=Syne:wght@700;800&display=swap');
  *{box-sizing:border-box;margin:0;padding:0}
  html,body,#root{background:#09090e;width:100%;height:100%;overflow:hidden}
  ::-webkit-scrollbar{width:6px}::-webkit-scrollbar-track{background:#0d0d0d}::-webkit-scrollbar-thumb{background:#2a4a32}
  @keyframes pulse{0%,100%{opacity:1}50%{opacity:.25}}
  @keyframes si{from{opacity:0;transform:translateY(5px)}to{opacity:1;transform:translateY(0)}}
  @keyframes dots{0%,100%{content:''}33%{content:'.'}66%{content:'..'}99%{content:'...'}}
  .thinking-indicator::after{content:'';animation:dots 1.2s steps(1) infinite}
  @keyframes spin{from{transform:rotate(0deg)}to{transform:rotate(360deg)}}
  .spinner{width:28px;height:28px;border:2px solid #1a3a22;border-top-color:#3dcc7a;border-radius:50%;animation:spin 0.9s linear infinite}
  .ent{animation:si .2s ease}
  .blink{animation:pulse 1.5s infinite}
  .btn{background:transparent;border:1px solid #2a5a3a;color:#3dcc7a;font-family:inherit;font-size:14px;padding:9px 16px;cursor:pointer;letter-spacing:2px;text-transform:uppercase;transition:all .15s;width:100%}
  .btn:hover{background:#0c1f16;border-color:#3dcc7a}.btn:disabled{opacity:.3;cursor:not-allowed}
  .bsm{padding:5px 12px;width:auto;font-size:13px}.bdng{border-color:#4a2222;color:#cc5555}.bdng:hover{background:#1a0a0a;border-color:#cc5555}
  .slot{border:1px solid #1e2e20;padding:14px;margin-bottom:10px;background:#0b0c0f}
  .slot.ok{border-color:#223a28}
  .uz{border:1px dashed #1e2e20;padding:12px;text-align:center;cursor:pointer;transition:all .15s;background:#0c0c10;display:block;margin-bottom:6px;font-size:14px;color:#3a6a4a}
  .uz:hover,.uz.ok{border-color:#3dcc7a;background:#0b160f;color:#3dcc7a}
  .tag{display:inline-block;padding:3px 8px;border-radius:2px;font-size:12px;font-weight:700;letter-spacing:1px}
  input[type=text],select{background:#0b0c0f;border:1px solid #1e2e20;color:#ccc;padding:6px 9px;font-size:14px;font-family:inherit;width:100%}
  input[type=number]{background:#0b0c0f;border:1px solid #1e2e20;color:#3dcc7a;padding:6px 9px;font-size:14px;font-family:inherit;width:100%}
  .sec{font-size:12px;color:#3a7a4a;letter-spacing:2px;margin:16px 0 9px;font-weight:600}
`;

const VERDICT_STYLE = {
  confirmed: { color: "#3dcc7a", icon: "✓" },
  rejected:  { color: "#cc4444", icon: "✗" },
  uncertain: { color: "#ccaa44", icon: "?" },
  pending:   { color: "#446688", icon: "○" },
};

export default function App() {
  const [slots, setSlots] = useState([
    { id: 0, exprFile: null, metaFile: null, name: "Dataset 1" },
    { id: 1, exprFile: null, metaFile: null, name: "Dataset 2" },
  ]);
  const [loaded,     setLoaded]     = useState([]);   // [{id, name, gene_count, sample_count, group_cols, group_col, groups}]
  const [groupMap,   setGroupMap]   = useState({});   // {dataset_id → chosen group_col}
  const [phase,      setPhase]      = useState("upload");
  const [log,        setLog]        = useState([]);
  const [hypotheses, setHypotheses] = useState([]);
  const [step,          setStep]          = useState(0);
  const [maxSteps,      setMaxSteps]      = useState(15);
  const [agentMode,     setAgentMode]     = useState("free");
  const [currentStatus, setCurrentStatus] = useState("");
  const [streamingText, setStreamingText] = useState("");
  const logEnd   = useRef(null);
  const abortRef = useRef(null);

  useEffect(() => { logEnd.current?.scrollIntoView({ behavior: "smooth" }); }, [log]);
  const addLog = useCallback(e => setLog(prev => [...prev, { ...e, id: Date.now() + Math.random() }]), []);

  const addSlot    = () => setSlots(p => [...p, { id: Date.now(), exprFile: null, metaFile: null, name: `Dataset ${p.length + 1}` }]);
  const removeSlot = id => setSlots(p => p.filter(s => s.id !== id));
  const updSlot    = (id, k, v) => setSlots(p => p.map(s => s.id === id ? { ...s, [k]: v } : s));

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
  };

  const runAgent = async () => {
    if (!loaded.length) return;
    setPhase("running"); setLog([]); setStep(0); setHypotheses([]); setCurrentStatus("Running pre-analysis..."); setStreamingText("");

    const controller = new AbortController();
    abortRef.current = controller;
    let currentStep = 0;

    try {
      const res = await fetch("http://localhost:8000/api/run", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ dataset_ids: loaded.map(d => d.id), group_cols: groupMap, max_steps: maxSteps, mode: agentMode }),
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
            if (entry.type === "thinking")        { setStep(++currentStep); setCurrentStatus(entry.text); setStreamingText(""); addLog(entry); continue; }
            if (entry.type === "thought_stream")  { flushSync(() => setStreamingText(prev => prev + entry.delta)); continue; }
            if (entry.type === "thought")         { flushSync(() => setStreamingText("")); }
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

  return (
    <div style={{ minHeight: "100vh", background: "#09090e", fontFamily: "'JetBrains Mono','Fira Code',monospace", color: "#ddd" }}>
      <style>{STYLES}</style>

      <div style={{ borderBottom: "1px solid #1a2e1a", padding: "16px 28px", display: "flex", alignItems: "center", gap: 12 }}>
        <div style={{ width: 8, height: 8, borderRadius: "50%", background: phase === "running" ? "#3dcc7a" : "#223322", boxShadow: phase === "running" ? "0 0 10px #3dcc7a" : "none" }} className={phase === "running" ? "blink" : ""} />
        <span style={{ fontFamily: "'Syne',sans-serif", fontSize: 20, fontWeight: 800, letterSpacing: 3, color: "#3dcc7a" }}>TRANSCRIPTOMIC AGENT</span>
        <span style={{ fontSize: 13, color: "#2a5a3a", letterSpacing: 2 }}>/ MULTI-DATASET · CROSS-COHORT</span>
        {phase === "running" && !currentStatus && <span style={{ marginLeft: "auto", fontSize: 14, color: "#3a7a4a" }}>STEP {step}/{maxSteps}</span>}
      </div>

      <div style={{ display: "flex", height: "calc(100vh - 58px)" }}>

        {/* LEFT PANEL */}
        <div style={{ width: 290, borderRight: "1px solid #1a2e1a", padding: "16px 14px", overflowY: "auto", flexShrink: 0 }}>
          <div className="sec">// DATASETS</div>

          {slots.map(slot => (
            <DatasetSlot key={slot.id} slot={slot} canRemove={slots.length > 1}
              onUpdate={(k, v) => updSlot(slot.id, k, v)} onRemove={() => removeSlot(slot.id)} />
          ))}

          <button className="btn bsm" style={{ marginBottom: 8, width: "100%" }} onClick={addSlot}>+ ADD DATASET</button>
          <button className="btn" style={{ marginBottom: 12 }} onClick={loadAll} disabled={!slots.some(s => s.exprFile && s.metaFile)}>
            LOAD DATA
          </button>

          {loaded.length > 0 && <>
            <div className="sec">// GROUP COLUMNS</div>
            {loaded.map(ds => (
              <div key={ds.id} style={{ marginBottom: 14 }}>
                <div style={{ fontSize: 14, color: "#4a9a6a", marginBottom: 5, fontWeight: 600 }}>{ds.name}</div>
                <select value={groupMap[ds.id] || ds.group_col} onChange={e => setGroupMap(prev => ({ ...prev, [ds.id]: e.target.value }))}>
                  {ds.group_cols.map(c => <option key={c} value={c}>{c}</option>)}
                </select>
                <div style={{ fontSize: 13, color: "#3a6a4a", marginTop: 5, lineHeight: 2 }}>
                  {ds.groups.map(g => <div key={g}>▸ {g}</div>)}
                </div>
                <div style={{ fontSize: 13, color: "#2a5a3a", marginTop: 4 }}>{ds.gene_count} genes · {ds.sample_count} samples</div>
              </div>
            ))}

            <div className="sec">// MODE</div>
            <div style={{ display: "flex", gap: 6, marginBottom: 12 }}>
              {["free", "hybrid"].map(m => (
                <button key={m} className="btn bsm" onClick={() => setAgentMode(m)}
                  style={{ flex: 1, borderColor: agentMode === m ? "#3dcc7a" : "#2a5a3a", color: agentMode === m ? "#3dcc7a" : "#2a5a3a", background: agentMode === m ? "#0b160f" : "transparent" }}>
                  {m.toUpperCase()}
                </button>
              ))}
            </div>

            <div className="sec">// MAX STEPS</div>
            <input type="number" value={maxSteps} min={3} max={30}
              onChange={e => setMaxSteps(parseInt(e.target.value))}
              style={{ marginBottom: 12 }} />

            <button className="btn" style={{ background: phase === "running" ? "#080e0a" : "transparent" }}
              onClick={phase === "running" ? () => abortRef.current?.abort() : runAgent}>
              {phase === "running" ? "■ STOP" : "▶ START AGENT"}
            </button>
          </>}
        </div>

        {/* LOG PANEL */}
        <div style={{ flex: 1, display: "flex", flexDirection: "column", overflow: "hidden" }}>

          {/* sticky status bar */}
          {currentStatus && (
            <div style={{ flexShrink: 0, display: "flex", alignItems: "center", gap: 10, padding: "10px 28px", borderBottom: "1px solid #1a2e1a", background: "#0b0f0b" }}>
              <div className="spinner" />
              <span className="thinking-indicator" style={{ fontSize: 13, color: "#3dcc7a", letterSpacing: 1 }}>{currentStatus}</span>
              {step > 0 && <span style={{ marginLeft: "auto", fontSize: 12, color: "#2a5a3a" }}>STEP {step}/{maxSteps}</span>}
            </div>
          )}

          <div style={{ flex: 1, overflowY: "auto", padding: "22px 28px" }}>
            {log.length === 0 && !currentStatus && (
              <div style={{ textAlign: "center", marginTop: 100 }}>
                <div style={{ fontSize: 36, opacity: .1, marginBottom: 14 }}>⬡</div>
                <div style={{ fontSize: 16, color: "#2a5a3a" }}>Load datasets and start the agent</div>
                <div style={{ fontSize: 13, color: "#1a3a22", marginTop: 8 }}>Backend: <code style={{color:"#2a5a3a"}}>uvicorn backend.main:app --reload</code></div>
              </div>
            )}
            {log.length === 0 && currentStatus && (
              <div style={{ textAlign: "center", marginTop: 120 }}>
                <div style={{ display: "flex", justifyContent: "center", marginBottom: 24 }}>
                  <div style={{ width: 48, height: 48, border: "2px solid #1a3a22", borderTopColor: "#3dcc7a", borderRadius: "50%", animation: "spin 0.9s linear infinite" }} />
                </div>
                <div style={{ fontSize: 15, color: "#3a7a4a", letterSpacing: 1 }}>{currentStatus}</div>
                <div style={{ fontSize: 12, color: "#1a3a22", marginTop: 8 }}>This may take a few seconds...</div>
              </div>
            )}
            {log.map(e => <LogEntry key={e.id} entry={e} />)}
            {streamingText && (
              <div className="ent" style={{ marginBottom: 12, borderLeft: "3px solid #224a2a", paddingLeft: 14 }}>
                <div style={{ display: "flex", alignItems: "flex-start", gap: 9 }}>
                  <span style={{ color: "#224a2a", fontSize: 15, marginTop: 2, flexShrink: 0 }}>◆</span>
                  <div style={{ fontSize: 14, color: "#5aaa7a", lineHeight: 1.7 }}>
                    {streamingText}<span className="blink">▋</span>
                  </div>
                </div>
              </div>
            )}
            <div ref={logEnd} />
          </div>
        </div>

        {/* HYPOTHESIS PANEL */}
        {(phase === "running" || hypotheses.length > 0) && (
          <div style={{ width: 270, borderLeft: "1px solid #1a2e1a", padding: "16px 14px", overflowY: "auto", flexShrink: 0, background: "#09090e" }}>
            <div className="sec">// HYPOTHESES</div>
            {hypotheses.length === 0 && <div style={{ fontSize: 13, color: "#2a4a2a" }}>Agent is formulating hypotheses...</div>}
            {hypotheses.map(h => {
              const vs = VERDICT_STYLE[h.status] || VERDICT_STYLE.pending;
              return (
                <div key={h.id} style={{ marginBottom: 12, padding: "10px 11px", background: "#0b0c0f", border: `1px solid ${vs.color}33`, borderRadius: 3 }}>
                  <div style={{ display: "flex", alignItems: "center", gap: 6, marginBottom: 6 }}>
                    <span className="tag" style={{ background: `${vs.color}22`, color: vs.color, border: `1px solid ${vs.color}44` }}>{h.id}</span>
                    <span style={{ fontSize: 14, color: vs.color }}>{vs.icon}</span>
                    <span style={{ fontSize: 12, color: vs.color, letterSpacing: 1, opacity: 0.8 }}>{h.status.toUpperCase()}</span>
                  </div>
                  <div style={{ fontSize: 13, color: "#7799bb", lineHeight: 1.7 }}>{h.text}</div>
                  {h.evidence.length > 0 && (
                    <div style={{ marginTop: 7, borderTop: "1px solid #1a2a3a", paddingTop: 7 }}>
                      {h.evidence.map((ev, i) => (
                        <div key={i} style={{ fontSize: 12, color: "#4a6677", lineHeight: 1.6, marginBottom: 4 }}>
                          <span style={{ color: "#3a5a6a" }}>krok {ev.step} [{ev.action}]</span> {ev.reasoning}
                          {ev.key_stats && Object.keys(ev.key_stats).length > 0 && (
                            <div style={{ marginTop: 2, paddingLeft: 8, borderLeft: "2px solid #1a3a2a" }}>
                              {Object.entries(ev.key_stats).map(([gene, s]) => (
                                <span key={gene} style={{ display: "inline-block", marginRight: 10, color: "#3a6a4a", fontSize: 11 }}>
                                  <b style={{ color: "#4a8a5a" }}>{gene}</b>{": "}
                                  {Object.entries(s).filter(([, v]) => v != null).map(([k, v]) =>
                                    `${k}=${typeof v === "number" ? (Math.abs(v) < 0.001 ? v.toExponential(2) : v.toPrecision(3)) : Array.isArray(v) ? v.join(",") : v}`
                                  ).join("  ")}
                                </span>
                              ))}
                            </div>
                          )}
                        </div>
                      ))}
                    </div>
                  )}
                </div>
              );
            })}
          </div>
        )}
      </div>
    </div>
  );
}
