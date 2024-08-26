import React from "react";

function SNPSummaryTable({ snps }) {
  // Provide a fallback for snps in case it's undefined or empty
  const snpsList = snps || [];

  return (
    <div>
      <h2>SNP Summary</h2>
      {snpsList.length > 0 ? (
        <table className="table table-striped table-bordered">
          <thead>
            <tr>
              <th>Position</th>
              <th>Reference Base</th>
              <th>Variant Base</th>
              <th>Annotation</th>
            </tr>
          </thead>
          <tbody>
            {snpsList.map((snp, index) => (
              <tr key={index}>
                <td>{snp.position}</td>
                <td>{snp.ref}</td>
                <td>{snp.qual}</td>
                <td>{snp.annotation || "N/A"}</td>{" "}
                {/* Default to "N/A" if annotation is not available */}
              </tr>
            ))}
          </tbody>
        </table>
      ) : (
        <p>No SNP data available.</p>
      )}
    </div>
  );
}

export default SNPSummaryTable;
