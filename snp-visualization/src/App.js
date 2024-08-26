import React, { useState, useEffect } from "react";
import GenomeViewer from "./Components/GenomeViewer";
import ZoomControls from "./Components/ZoomControls";
import VarianceBoxPlot from "./Components/VarianceBoxPlot";
import FullPageGraph from "./Components/FullPageGraph";
import SequenceVisualizer from "./Components/SequenceVisualizer";
import SNPSummaryTable from "./Components/SNPSummaryTable";
import "bootstrap/dist/css/bootstrap.min.css";

function App() {
  const [snps, setSnps] = useState([]);
  const [zoomLevel, setZoomLevel] = useState(1);
  const [sequence, setSequence] = useState("ACTGGTACGTAGCTAGGCTAGCTAGCTAGC"); // Sample sequence data

  useEffect(() => {
    // Fetch the SNP data from the Flask backend
    fetch("http://127.0.0.1:8000/get-snps?vcf_file=output_variants.vcf")
      .then((response) => response.json())
      .then((data) => setSnps(data))
      .catch((error) => console.error("Error fetching SNP data:", error));
  }, []);

  return (
    <div
      style={{
        width: "100%",
        height: "100%",
        display: "flex",
        flexDirection: "column",
      }}
    >
      <header style={{ textAlign: "center", padding: "20px" }}>
        <h1>Genome Viewer</h1>
      </header>

      <main style={{ flex: 1, display: "flex", flexDirection: "row" }}>
        {/* Left Sidebar for Zoom Controls */}
        <aside style={{ width: "20%", padding: "10px" }}>
          <ZoomControls zoomLevel={zoomLevel} setZoomLevel={setZoomLevel} />
          <VarianceBoxPlot data={snps} />
        </aside>

        {/* Main Content Area for Genome Viewer, Sequence Visualizer, and FullPageGraph */}
        <section style={{ width: "80%", padding: "10px" }}>
          <GenomeViewer snps={snps} zoomLevel={zoomLevel} />
          <div style={{ marginTop: "20px" }}>
            <SequenceVisualizer sequence={sequence} snps={snps} />
          </div>
          <div style={{ marginTop: "20px" }}>
            <FullPageGraph data={snps} />
          </div>
          <SNPSummaryTable snps={snps} />
        </section>
      </main>
    </div>
  );
}

export default App;
