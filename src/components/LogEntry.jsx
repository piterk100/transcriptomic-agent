import { useState } from "react";

const ACTION_COLORS = {
  dataset_summary:           "#2a8a5a",
  top_variable_genes:        "#2a6a8a",
  differential_expression:   "#8a3a5a",
  gene_expression_by_group:  "#4a7a4a",
  nonlinear_rule:            "#8a5a2a",
  contextual_modules:        "#7a6a2a",
  pathway_enrichment:        "#5a8a3a",
  batch_detection:           "#7a5a8a",
  subgroup_discovery:        "#3a7a8a",
  gene_network_hub:          "#8a8a3a",
  cross_dataset_de:          "#3dff8a",
  cross_dataset_correlation: "#3dddff",
  invariant_axis:            "#ffdd3d",
  cross_dataset_rewiring:    "#ff8a3d",
  execute_code:              "#cc88ff",
};

const VERDICT_STYLE = {
  confirmed: { color: "#3dcc7a", icon: "✓", label: "CONFIRMED" },
  rejected:  { color: "#cc4444", icon: "✗", label: "REJECTED"  },
  uncertain: { color: "#ccaa44", icon: "?", label: "UNCERTAIN"  },
  pending:   { color: "#5588aa", icon: "○", label: "PENDING"    },
};

const TYPE_STYLES = {
  thinking:           { c: "#1e3a22", icon: "○" },
  thought:            { c: "#224a2a", icon: "◆" },
  result:             { c: "#1a2e1e", icon: "▸" },
  error:              { c: "#3a1818", icon: "✗" },
  done:               { c: "#224a2a", icon: "✓" },
  code:               { c: "#221a3a", icon: "⌥" },
  hypothesis_propose: { c: "#1a2e4a", icon: "H" },
  hypothesis_eval:    { c: "#221a3a", icon: "H" },
  seed:               { c: "#1a2a4a", icon: "⬡" },
};

export default function LogEntry({ entry }) {
  const [expanded, setExpanded] = useState(false);
  const ts = TYPE_STYLES[entry.type] || { c: "#1e3a22", icon: "·" };
  const ac = ACTION_COLORS[entry.action] || "#4a7a5a";

  return (
    <div className="ent" style={{ marginBottom: 12, borderLeft: `3px solid ${ts.c}`, paddingLeft: 14 }}>
      <div
        style={{ display: "flex", alignItems: "flex-start", gap: 9, cursor: entry.type === "result" ? "pointer" : "default" }}
        onClick={() => entry.type === "result" && setExpanded(e => !e)}
      >
        <span style={{ color: ts.c, fontSize: 13, marginTop: 2, flexShrink: 0 }}>{ts.icon}</span>
        <div style={{ flex: 1, minWidth: 0 }}>

          {entry.type === "thinking" && (
            <span style={{ fontSize: 12, color: "#2a5a3a" }}>{entry.text}</span>
          )}

          {entry.type === "thought" && (
            <div style={{ fontSize: 13, color: "#5aaa7a", lineHeight: 1.7 }}>{entry.text}</div>
          )}

          {entry.type === "code" && (
            <div>
              <div style={{ fontSize: 11, color: "#8866cc", marginBottom: 5 }}>⌥ agent writing custom code</div>
              <pre style={{ padding: 10, background: "#0e0e1a", border: "1px solid #2a2a4a", fontSize: 11, color: "#aa88ee", overflowX: "auto", maxHeight: 150, lineHeight: 1.6, borderRadius: 3 }}>
                {entry.code}
              </pre>
            </div>
          )}

          {entry.type === "result" && (
            <div>
              <div style={{ display: "flex", alignItems: "center", gap: 7, flexWrap: "wrap" }}>
                <span
                  className="tag"
                  style={{ background: `${ac}22`, color: ac, border: `1px solid ${ac}44`, fontWeight: entry.isCross || entry.isDynamic ? 700 : 500, fontSize: 11 }}
                >
                  {entry.isCross ? "⬡ " : entry.isDynamic ? "⌥ " : ""}{entry.action}
                </span>
                <span style={{ fontSize: 12, color: "#5a9a6a", flex: 1 }}>{entry.summary}</span>
                <span style={{ fontSize: 11, color: "#2a5a3a" }}>{expanded ? "▲" : "▼"}</span>
              </div>
              {expanded && (
                <pre style={{ marginTop: 9, padding: 11, background: "#0a120a", border: "1px solid #1a2e1a", fontSize: 11, color: "#3a8a4a", overflowX: "auto", maxHeight: 300, overflowY: "auto", lineHeight: 1.6, borderRadius: 3 }}>
                  {JSON.stringify(entry.result, null, 2).slice(0, 4000)}
                </pre>
              )}
            </div>
          )}

          {entry.type === "hypothesis_propose" && (
            <div style={{ padding: "8px 10px", background: "#0d1a2a", border: "1px solid #2a3a5a", borderRadius: 3 }}>
              <div style={{ display: "flex", alignItems: "center", gap: 8, marginBottom: 5 }}>
                <span className="tag" style={{ background: "#2a4a7a", color: "#88aadd", fontSize: 10, border: "none" }}>
                  {entry.hypothesis.id}
                </span>
                <span style={{ fontSize: 11, color: "#5588aa", letterSpacing: 1 }}>NEW HYPOTHESIS</span>
              </div>
              <div style={{ fontSize: 13, color: "#88aabb", lineHeight: 1.7 }}>{entry.hypothesis.text}</div>
            </div>
          )}

          {entry.type === "hypothesis_eval" && (() => {
            const vs = VERDICT_STYLE[entry.hypothesis.status] || VERDICT_STYLE.uncertain;
            return (
              <div style={{ padding: "8px 10px", background: "#0d0f0d", border: `1px solid ${vs.color}33`, borderRadius: 3 }}>
                <div style={{ display: "flex", alignItems: "center", gap: 8, marginBottom: 5 }}>
                  <span className="tag" style={{ background: `${vs.color}22`, color: vs.color, fontSize: 10, border: `1px solid ${vs.color}44` }}>
                    {entry.hypothesis.id}
                  </span>
                  <span style={{ fontSize: 12, color: vs.color, letterSpacing: 1 }}>{vs.icon} {vs.label}</span>
                </div>
                {entry.reasoning && (
                  <div style={{ fontSize: 12, color: "#6688aa", lineHeight: 1.7 }}>{entry.reasoning}</div>
                )}
              </div>
            );
          })()}

          {entry.type === "seed" && (
            <div style={{ padding: "8px 10px", background: "#0d0f18", border: "1px solid #2a3a5a", borderRadius: 3 }}>
              <div style={{ fontSize: 11, color: "#5577aa", letterSpacing: 1, marginBottom: 5 }}>⬡ STATISTICAL PRE-ANALYSIS</div>
              {entry.summary
                ? <div style={{ fontSize: 11, color: "#4466aa", lineHeight: 1.7, whiteSpace: "pre-line" }}>{entry.summary}</div>
                : <div style={{ fontSize: 11, color: "#334466" }}>No common groups across datasets — agent starting without seed hypotheses.</div>
              }
            </div>
          )}

          {entry.type === "error" && (
            <span style={{ fontSize: 13, color: "#cc5555" }}>{entry.text}</span>
          )}

          {entry.type === "done" && (
            <span style={{ fontSize: 14, color: "#3dcc7a", fontWeight: 600 }}>{entry.text}</span>
          )}
        </div>
      </div>
    </div>
  );
}
