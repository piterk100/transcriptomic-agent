export async function setGroupMappings(mappings) {
  const response = await fetch("/api/group_mappings", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ mappings }),
  });
  return response.json();
}

export async function uploadDegDataset(file, name, groupA, groupB) {
  const formData = new FormData();
  formData.append("file", file);
  formData.append("name", name);
  formData.append("groupA", groupA);
  formData.append("groupB", groupB);
  const response = await fetch("/api/datasets/upload_deg", {
    method: "POST",
    body: formData,
  });
  return response.json();
}
