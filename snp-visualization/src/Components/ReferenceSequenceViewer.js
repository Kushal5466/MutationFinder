import React from "react";

function ReferenceSequenceViewer({ sequence, annotations, zoomLevel }) {
  // You can use zoomLevel to adjust the displayed sequence length
  const displayedSequence = sequence.slice(0, sequence.length / zoomLevel);

  return (
    <div>
      <h2>Reference Sequence</h2>
      <p>{displayedSequence}</p>
      <h3>Annotations</h3>
      <ul>
        {annotations.map((annotation, index) => (
          <li key={index}>
            {annotation.name}: {annotation.position}
          </li>
        ))}
      </ul>
    </div>
  );
}

export default ReferenceSequenceViewer;
