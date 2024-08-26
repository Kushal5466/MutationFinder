import React from "react";

// Functional component to display a summary table of SNPs
function SNPSummaryTable({ snps }) {
  // Provide a fallback for snps in case it's undefined or empty
  const snpsList = snps || [];

  return (
    <div>
      {/* Section header for the SNP summary */}
      <h2>SNP Summary</h2>

      {/* Conditional rendering: if SNP data is available, display the table; otherwise, show a message */}
      {snpsList.length > 0 ? (
        <table className="table table-striped table-bordered">
          <thead>
            <tr>
              {/* Table headers to describe the columns of the SNP data */}
              <th>Position</th>
              <th>Reference Base</th>
              <th>Variant Base</th>
              <th>Annotation</th>
            </tr>
          </thead>
          <tbody>
            {/* Iterate over the SNP list and render each SNP in its own row */}
            {snpsList.map((snp, index) => (
              <tr key={index}>
                {/* Display the SNP data in each cell */}
                <td>{snp.position}</td>
                <td>{snp.ref}</td>
                <td>{snp.qual}</td>
                {/* Show "N/A" if the annotation is not available */}
                <td>{snp.annotations || "N/A"}</td>
              </tr>
            ))}
          </tbody>
        </table>
      ) : (
        // Fallback message if no SNP data is provided
        <p>No SNP data available.</p>
      )}
    </div>
  );
}

// Export the component to be used in other parts of the application
export default SNPSummaryTable;
