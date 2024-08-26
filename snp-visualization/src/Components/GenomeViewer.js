import React, { useRef, useEffect } from "react";
import * as d3 from "d3";

const GenomeViewer = ({ snps, zoomLevel }) => {
  const svgRef = useRef();
  const width = 1000; // SVG width
  const height = 100; // SVG height

  useEffect(() => {
    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove(); // Clear previous renders

    // Scale for positioning SNPs based on zoom level
    const xScale = d3
      .scaleLinear()
      .domain([0, 1000000 / zoomLevel]) // Adjust domain based on zoom
      .range([0, width]);

    // Draw genome line
    svg
      .append("line")
      .attr("x1", 0)
      .attr("y1", height / 2)
      .attr("x2", width)
      .attr("y2", height / 2)
      .attr("stroke", "black")
      .attr("stroke-width", 2);

    // Plot SNPs as circles on the genome line
    snps.forEach((snp) => {
      svg
        .append("circle")
        .attr("cx", xScale(snp.position))
        .attr("cy", height / 2)
        .attr("r", 5)
        .attr("fill", "red")
        .append("title")
        .text(`Position: ${snp.position}, Mutation: ${snp.ref} > ${snp.alt}`);
    });
  }, [snps, zoomLevel]);

  return <svg ref={svgRef} width={width} height={height}></svg>;
};

export default GenomeViewer;
