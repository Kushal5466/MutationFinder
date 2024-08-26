import React, { useEffect, useState } from "react";
import Plot from "react-plotly.js";

function FullPagePlotlyGraph({ data }) {
  const positions = data.map((entry) => entry.position);
  const qualValues = data.map((entry) => parseFloat(entry.qual));
  const chromosome = data.map((entry) => entry.chromosome);
  const labels = data.map(
    (entry) => `QUAL: ${entry.chromosome}, Position: ${entry.position}`
  );
  const [windowSize, setWindowSize] = useState({
    width: window.innerWidth,
    height: window.innerHeight,
  });

  useEffect(() => {
    const handleResize = () => {
      setWindowSize({
        width: window.innerWidth,
        height: window.innerHeight,
      });
    };

    window.addEventListener("resize", handleResize);

    // Cleanup event listener on component unmount
    return () => window.removeEventListener("resize", handleResize);
  }, []);

  return (
    <Plot
      data={[
        {
          x: positions,
          y: qualValues,
          type: "scatter",
          mode: "lines+markers",
          marker: { color: "red" },
        },
        // { type: "bar", x: positions, y: qualValues },
      ]}
      layout={{
        width: windowSize.width,
        height: windowSize.height,
        title: "Full-Page Plotly Graph",
      }}
      style={{ width: "100vw", height: "100vh" }}
      config={{ responsive: true }}
    />
  );
}

export default FullPagePlotlyGraph;
