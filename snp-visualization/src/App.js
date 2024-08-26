import React, { useState, useEffect } from "react";
import { BrowserRouter as Router, Route, Routes, Link } from "react-router-dom";
import GenomeViewer from "./Components/GenomeViewer";
import VarianceBoxPlot from "./Components/VarianceBoxPlot";
import SNPSummaryTable from "./Components/SNPSummaryTable";
import ZoomControls from "./Components/ZoomControls"; // Import ZoomControls component for zoom functionality
import "bootstrap/dist/css/bootstrap.min.css"; // Import Bootstrap CSS for styling

function App() {
  // State to hold SNP data fetched from the backend
  const [snps, setSnps] = useState([]);
  // State to manage the zoom level for the GenomeViewer component
  const [zoomLevel, setZoomLevel] = useState(1);

  // useEffect hook to fetch SNP data from the backend API when the component mounts
  useEffect(() => {
    fetch("http://127.0.0.1:8000/get-snps?vcf_file=output_variants.vcf")
      .then((response) => response.json()) // Parse the response as JSON
      .then((data) => setSnps(data)) // Update the snps state with the fetched data
      .catch((error) => console.error("Error fetching SNP data:", error)); // Log any errors to the console
  }, []); // Empty dependency array means this effect runs once on component mount

  return (
    <Router>
      <div
        className="container-fluid"
        style={{ height: "100vh", display: "flex", flexDirection: "column" }}
      >
        {/* Navbar Section */}
        <header style={{ textAlign: "center", padding: "20px" }}>
          <h1>Genome Viewer App</h1>
          <nav>
            {/* Links to navigate between different routes */}
            <Link to="/genome-viewer" style={{ margin: "0 10px" }}>
              Genome Viewer
            </Link>
            <Link to="/variance-box-plot" style={{ margin: "0 10px" }}>
              Variance Box Plot
            </Link>
            <Link to="/snp-summary-table" style={{ margin: "0 10px" }}>
              SNP Summary Table
            </Link>
          </nav>
        </header>

        <Routes>
          {/* Route for Genome Viewer with Zoom Controls */}
          <Route
            path="/genome-viewer"
            element={
              <div>
                {/* Integrate ZoomControls component to control zoom level */}
                <ZoomControls
                  zoomLevel={zoomLevel}
                  setZoomLevel={setZoomLevel}
                />
                {/* Pass snps and zoomLevel as props to GenomeViewer component */}
                <GenomeViewer
                  snps={snps}
                  zoomLevel={zoomLevel}
                  setZoomLevel={setZoomLevel}
                />
              </div>
            }
          />

          {/* Route for Variance Box Plot */}
          <Route
            path="/variance-box-plot"
            element={<VarianceBoxPlot data={snps} />}
          />

          {/* Route for SNP Summary Table */}
          <Route
            path="/snp-summary-table"
            element={<SNPSummaryTable snps={snps} />}
          />

          {/* Default route (/) also showing Genome Viewer with Zoom Controls */}
          <Route
            path="/"
            element={
              <div>
                {/* Integrate ZoomControls component to control zoom level */}
                <ZoomControls
                  zoomLevel={zoomLevel}
                  setZoomLevel={setZoomLevel}
                />
                {/* Pass snps and zoomLevel as props to GenomeViewer component */}
                <GenomeViewer
                  snps={snps}
                  zoomLevel={zoomLevel}
                  setZoomLevel={setZoomLevel}
                />
              </div>
            }
          />
        </Routes>
      </div>
    </Router>
  );
}

export default App;
