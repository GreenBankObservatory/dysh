// download.js
window.addEventListener("load", () => {
  const links = document.querySelectorAll(
    '.dropdown-download-buttons a[href$=".ipynb"]'
  );

  for (const link of links) {
    link.setAttribute("download", "");
  }
});
