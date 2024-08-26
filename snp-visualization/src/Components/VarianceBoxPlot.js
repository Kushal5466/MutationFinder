// src/VarianceScatterPlot.js
import React from "react";
import Plot from "react-plotly.js";

const VarianceScatterPlot = ({ data }) => {
  const positions = data.map((entry) => entry.position);
  const qualValues = data.map((entry) => parseFloat(entry.qual));
  const labels = data.map(
    (entry) => `QUAL: ${entry.qual}, Position: ${entry.position}`
  );

  return (
    <Plot
      data={[
        {
          x: positions,
          y: qualValues,
          mode: "lines+markers",
          marker: { color: "blue" },
          line: { shape: "linear" },
          text: labels,
          hoverinfo: "text", // Displays the QUAL value and position on hover
        },
      ]}
      layout={{
        title: "QUAL Values Scatter Plot",
        xaxis: {
          title: "Position",
        },
        yaxis: {
          title: "QUAL Value",
        },
        hovermode: "closest",
        dragmode: "zoom", // Enables zooming
      }}
      config={{
        responsive: true,
        scrollZoom: true, // Enables zooming with the mouse wheel
      }}
    />
  );
};

export default VarianceScatterPlot;
