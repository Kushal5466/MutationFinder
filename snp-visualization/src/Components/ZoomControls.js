import React from "react";

const ZoomControls = ({ zoomLevel, setZoomLevel }) => {
  const handleZoomIn = () => {
    setZoomLevel((prevZoomLevel) => prevZoomLevel * 2); // Zoom in
  };

  const handleZoomOut = () => {
    setZoomLevel((prevZoomLevel) => Math.max(1, prevZoomLevel / 2)); // Zoom out, but don't go below 1x zoom
  };

  return (
    <div>
      <button onClick={handleZoomIn}>Zoom In</button>
      <button onClick={handleZoomOut}>Zoom Out</button>
    </div>
  );
};

export default ZoomControls;
