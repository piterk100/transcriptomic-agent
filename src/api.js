export async function setGroupMappings(mappings) {
  const response = await fetch("/api/group_mappings", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ mappings }),
  });
  return response.json();
}
