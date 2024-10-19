"use strict";

// This is not a full .obj parser.
// see http://paulbourke.net/dataformats/obj/

function parseOBJ(text) {
  // because indices are base 1 let's just fill in the 0th data
  const objPositions = [[0, 0, 0]];
  const objTexcoords = [[0, 0]];
  const objNormals = [[0, 0, 0]];

  // same order as `f` indices
  const objVertexData = [
    objPositions,
    objTexcoords,
    objNormals,
  ];

  // same order as `f` indices
  let webglVertexData = [
    [],   // positions
    [],   // texcoords
    [],   // normals
  ];

  function newGeometry() {
    // If there is an existing geometry and it's
    // not empty then start a new one.
    if (geometry && geometry.data.position.length) {
      geometry = undefined;
    }
    setGeometry();
  }

  function addVertex(vert) {
    const ptn = vert.split('/');
    ptn.forEach((objIndexStr, i) => {
      if (!objIndexStr) {
        return;
      }
      const objIndex = parseInt(objIndexStr);
      const index = objIndex + (objIndex >= 0 ? 0 : objVertexData[i].length);
      webglVertexData[i].push(...objVertexData[i][index]);
    });
  }

  const keywords = {
    v(parts) {
      objPositions.push(parts.map(parseFloat));
    },
    vn(parts) {
      objNormals.push(parts.map(parseFloat));
    },
    vt(parts) {
      // should check for missing v and extra w?
      objTexcoords.push(parts.map(parseFloat));
    },
    f(parts) {
      const numTriangles = parts.length - 2;
      for (let tri = 0; tri < numTriangles; ++tri) {
        addVertex(parts[0]);
        addVertex(parts[tri + 1]);
        addVertex(parts[tri + 2]);
      }
    },
  };

  const keywordRE = /(\w*)(?: )*(.*)/;
  const lines = text.split('\n');
  for (let lineNo = 0; lineNo < lines.length; ++lineNo) {
    const line = lines[lineNo].trim();
    if (line === '' || line.startsWith('#')) {
      continue;
    }
    const m = keywordRE.exec(line);
    if (!m) {
      continue;
    }
    const [, keyword, unparsedArgs] = m;
    const parts = line.split(/\s+/).slice(1);
    const handler = keywords[keyword];
    if (!handler) {
      console.warn('unhandled keyword:', keyword);  // eslint-disable-line no-console
      continue;
    }
    handler(parts, unparsedArgs);
  }

  return {
    position: webglVertexData[0],
    texcoord: webglVertexData[1],
    normal: webglVertexData[2],
  };
}

function getIndices(array)
{
    let ret = []
    for (let i = 0; i < array.length; i++){
        ret.push(i)
    } 
    return ret
}

async function main() {
  // Get A WebGL context
  /** @type {HTMLCanvasElement} */
    const canvas = document.getElementById("canvas");

  const gl = canvas.getContext("webgl");
  if (!gl) {
    return;
  }


  // compiles and links the shaders, looks up attribute and uniform locations
 // const meshProgramInfo = webglUtils.createProgramInfo(gl, [vs, fs]);
  const response = await fetch('./data/cells/stromal_like.obj');  
  const text = await response.text();
  const data = parseOBJ(text);

  // Because data is just named arrays like this
  //
  // {
  //   position: [...],
  //   texcoord: [...],
  //   normal: [...],
  // }
  //
  // and because those names match the attributes in our vertex
  // shader we can pass it directly into `createBufferInfoFromArrays`
  // from the article "less code more fun".

  // create a buffer for each array by calling
  // gl.createBuffer, gl.bindBuffer, gl.bufferData
  //const bufferInfo = webglUtils.createBufferInfoFromArrays(gl, data);
  
  // position buffer
  function createBuffer(gl, data) {
      const bin_buffer = gl.createBuffer()

      const bin_normal = new Float32Array(data.normal)
//      const tex_i = getIndices(data.normal)
      gl.bindBuffer(gl.ARRAY_BUFFER, bin_buffer)
      gl.bufferData(gl.ARRAY_BUFFER, bin_normal, gl.STATIC_DRAW)

      const tex_buffer = gl.createBuffer()
      const bin_texture = new Float32Array(data.texcoord)
//      const tex_i = getIndices(data.texcoord)
      gl.bindBuffer(gl.ARRAY_BUFFER, tex_buffer)
      gl.bufferData(gl.ARRAY_BUFFER, bin_texture, gl.STATIC_DRAW)
      // This could be gemini hallucinating, but i think we need to use typed
      // arrays so its actually binary data
      
      const pos_buffer = gl.createBuffer()
      const bin_positions = new Float32Array(data.position)
 //     const pos_i = getIndices(data.position)
      gl.bindBuffer(gl.ARRAY_BUFFER, pos_buffer)
      gl.bufferData(gl.ARRAY_BUFFER, bin_positions, gl.STATIC_DRAW)
  }
  createBuffer(gl, data);

  const cameraTarget = [0, 0, 0];
  const cameraPosition = [0, 0, 4];
  const zNear = 0.1;
  const zFar = 50;

  function degToRad(deg) {
    return deg * Math.PI / 180;
  }

  function render(time) {
    time *= 0.001;  // convert to seconds

//    webglUtils.resizeCanvasToDisplaySize(gl.canvas);

      // non-dynamic canvas size
    gl.canvas.width = window.innerWidth;
    gl.canvas.height = window.innerHeight;
    gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
    gl.enable(gl.DEPTH_TEST);
    gl.enable(gl.CULL_FACE);

           // Setup the pointer to our attribute data (the triangles)
      gl.enableVertexAttribArray(gl.positionLocation);
      gl.vertexAttribPointer(gl.positionLocation, 3, gl.FLOAT, false, 0, 0);

      // lets just make them red for now
      const red = new Float32Array([0.0, 0.0, 1.0, 1.0])
      // Setup the color uniform that will be shared across all triangles
      gl.uniform4fv(gl.colorLocation, red);

      // Draw the triangles to the screen
      gl.drawArrays(gl.TRIANGLES, 0, 6);


    requestAnimationFrame(render);
  }
  requestAnimationFrame(render);
}

main();
