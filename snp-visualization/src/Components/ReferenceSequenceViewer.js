import React from "react";

// Functional component to display a reference sequence and its annotations
function ReferenceSequenceViewer({ sequence, annotations, zoomLevel }) {
  // Used the zoomLevel prop to adjust the length of the sequence displayed to the user
  // The displayed sequence is a subset of the full sequence, determined by the zoomLevel
  const displayedSequence = sequence.slice(0, sequence.length / zoomLevel);

  return (
    <div>
      {/* Section header for the reference sequence */}
      <h2>Reference Sequence</h2>
      {/* Display the adjusted (zoomed) sequence */}
      <p>{displayedSequence}</p>

      {/* Section header for annotations */}
      <h3>Annotations</h3>
      {/* Display each annotation as a list item */}
      <ul>
        {annotations.map((annotation, index) => (
          // Key is provided to ensure each list item is uniquely identifiable by React
          <li key={index}>
            {/* Display the annotation's name and its corresponding position */}
            {annotation.name}: {annotation.position}
          </li>
        ))}
      </ul>
    </div>
  );
}

// Export the component for use in other parts of the application
export default ReferenceSequenceViewer;
