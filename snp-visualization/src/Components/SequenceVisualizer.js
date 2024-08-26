// src/components/SequenceVisualizer.js

import React, { useEffect, useRef } from "react";
import * as d3 from "d3";
import "./SequenceVisualizer.css";

const SequenceVisualizer = ({ sequence, snps }) => {
  const svgRef = useRef();
  console.log(snps);
  useEffect(() => {
    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove(); // Clear previous renders

    const width = 800;
    const height = 200;
    const margin = { top: 20, right: 30, bottom: 30, left: 40 };

    svg
      .attr("viewBox", `0 0 ${width} ${height}`)
      .style("background", "#f0f0f0");

    const xScale = d3
      .scaleLinear()
      .domain([1, sequence.length])
      .range([margin.left, width - margin.right]);

    const xAxis = d3
      .axisBottom(xScale)
      .ticks(20)
      .tickSize(-height + margin.top + margin.bottom);

    svg
      .append("g")
      .attr("transform", `translate(0,${height - margin.bottom})`)
      .call(xAxis);

    svg
      .append("g")
      .selectAll("rect")
      .data(sequence.split(""))
      .enter()
      .append("rect")
      .attr("x", (d, i) => xScale(i + 1))
      .attr("y", margin.top)
      .attr("width", (width - margin.left - margin.right) / sequence.length)
      .attr("height", height - margin.top - margin.bottom)
      .attr("fill", (d, i) => {
        const snp = snps.find((snp) => snp.position === i + 1);
        return snp ? "red" : "gray";
      });

    svg
      .append("g")
      .selectAll("text")
      .data(sequence.split(""))
      .enter()
      .append("text")
      .attr("x", (d, i) => xScale(i + 1))
      .attr("y", margin.top - 5)
      .attr("text-anchor", "middle")
      .attr("font-size", "10px")
      .text((d) => d);

    // Add SNP markers
    snps.forEach((snp) => {
      svg
        .append("rect")
        .attr("x", xScale(snp.position))
        .attr("y", margin.top)
        .attr("width", (width - margin.left - margin.right) / sequence.length)
        .attr("height", height - margin.top - margin.bottom)
        .attr("fill", "red")
        .attr("opacity", 0.5);
    });
  }, [sequence, snps]);

  return (
    <div className="sequence-visualizer">
      <h4>Reference Sequence</h4>
      <svg ref={svgRef}></svg>
    </div>
  );
};

export default SequenceVisualizer;
