var button = document.getElementById("downloadBtn");
button.onclick = function() {
    var imageOut = document.getElementById('svgoutput');

    domtoimage.toPng(imageOut, {
        width: document.getElementById('width').getAttribute("value"),
        height: document.getElementById('height').getAttribute("value")
      }).then(function (blob) {
        window.saveAs(blob, 'output.png');
    })
    .catch(function (error) {
        console.error('oops, something went wrong!', error);
    });
      
    
    
}