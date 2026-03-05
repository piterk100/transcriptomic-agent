export default function DatasetSlot({ slot, onUpdate, onRemove, canRemove }) {
  return (
    <div className={`slot ${slot.exprFile && slot.metaFile ? "ok" : ""}`}>
      <div style={{ display: "flex", gap: 5, marginBottom: 7 }}>
        <input
          type="text"
          value={slot.name}
          onChange={e => onUpdate("name", e.target.value)}
          style={{ flex: 1 }}
        />
        {canRemove && (
          <button className="btn bsm bdng" onClick={onRemove}>✕</button>
        )}
      </div>
      <label className={`uz ${slot.exprFile ? "ok" : ""}`}>
        <input type="file" accept=".csv" style={{ display: "none" }} onChange={e => onUpdate("exprFile", e.target.files[0])} />
        <span style={{ fontSize: 12, color: slot.exprFile ? "#3dcc7a" : "#3a6a4a" }}>
          {slot.exprFile ? `✓ ${slot.exprFile.name.slice(0, 24)}` : "+  macierz ekspresji"}
        </span>
      </label>
      <label className={`uz ${slot.metaFile ? "ok" : ""}`}>
        <input type="file" accept=".csv" style={{ display: "none" }} onChange={e => onUpdate("metaFile", e.target.files[0])} />
        <span style={{ fontSize: 12, color: slot.metaFile ? "#3dcc7a" : "#3a6a4a" }}>
          {slot.metaFile ? `✓ ${slot.metaFile.name.slice(0, 24)}` : "+  metadane (pheno)"}
        </span>
      </label>
    </div>
  );
}
