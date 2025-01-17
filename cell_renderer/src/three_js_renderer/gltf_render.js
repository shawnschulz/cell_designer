import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import * as THREE from 'three';

const scene = new THREE.Scene();
scene.background = new THREE.Color(0xf0f0f0);

// Camera setup
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.z = 5;

// Renderer setup
const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

// Lighting
const ambientLight = new THREE.AmbientLight(0xffffff, 0.5);
scene.add(ambientLight);

const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
directionalLight.position.set(0, 1, 0);
scene.add(directionalLight);

// Orbit Controls
const controls = new OrbitControls( camera, renderer.domElement );

//controls.update() must be called after any manual changes to the camera's transform
camera.position.set( 0, 1, 5);
controls.update();

// GLTF Loader
const loader = new GLTFLoader();
let model = new THREE.Object3D();

loader.load(
    './data/cells/stromal_like.glb', (gltf) => {
        model = gltf.scene;
        model.name = 'stromal_cell';
        scene.add(model);
        model.position.set(0,0,0);
    }
);

//loader.load(
//    // Resource URL - replace with your GLTF file path
//    './data/cells/stromal_like.glb',
//    
//    // Called when the resource is loaded
//    function (gltf) {
//        scene.add(gltf.scene);
//    },
//    
//    // Called while loading is progressing
//    function (xhr) {
//        console.log((xhr.loaded / xhr.total * 100) + '% loaded');
//    },
//    
//    // Called when loading has errors
//    function (error) {
//        console.error('An error happened', error);
//    }
//);

// Make 3 different colored cubes for the demo
//const geometry = new THREE.BoxGeometry( 1, 1, 1 ); 
//const material = new THREE.MeshBasicMaterial( {color: 0x00ff00} ); 
//const cube = new THREE.Mesh( geometry, material ); 


function placeOnSurface(gltf_object, protein_render_level=100, receptor_abundance_precent, receptor_name, receptor_color) {
    // Take a source object, get a random set of vertices, then place
    // the surface receptor gltf_object at the vertex
    // Can go with this for now, however if our render targets become more complex
    // may want to consider using threejs's InstancedMesh class instead of the default Mesh
    const geometry = new THREE.BoxGeometry( 1, 1, 1 ); 
    const material = new THREE.MeshBasicMaterial( {color: 0x00ff00} ); 
    const cube = new THREE.Mesh( geometry, material ); 
    scene.add( cube );
    cube.position.set(1.62,1.62,1.62);
    cube.scale.set(0.25,0.25,0.25);
    // How do we access the geomtery of this object?
    var object = scene.getObjectByName( "stromal_cell" );
    const positionAttribute = object.geometry.getAttribute('position');
    const vertex = new THREE.Vector3();
    // var verticies = gltf_object.geometry.position.array
    for ( let vertexIndex = 0; vertexIndex < positionAttribute.count; vertexIndex ++ ) {
        vertex.fromBufferAttribute( positionAttribute, vertexIndex );
    }
    // test by just placing something at the first vertex

    // Want the world position of the vertex
    gltf_object.localToWorld(vertex);
    cube.position = vertex;

    // Math.random()

}


// Rendering loop
function animate() {
    requestAnimationFrame(animate);

	// required if controls.enableDamping or controls.autoRotate are set to true
	controls.update();

    renderer.render(scene, camera);
}
animate();
const gltf_object = scene.getObjectByName( "stromal_cell" );
placeOnSurface(gltf_object);

function getProteins() {
    // Should do an API request to the backend to get the proteins
    // the user wants to add to the cell. This should create a mapping
    // of the protein name to a color and shape, we can save the mapping
    // to an SQL database and either in this or separate function 
    // go through all the mappings and create receptor objects and place
    // on surface
}

function getReceptorObject(color, shape) {
    // Generates the requested receptor mesh. 
    // Gonna start with reading the shape from a file, but
    // could change this to preload shapes
    // Color should either programatically create the material
    // or read the color material from a file
}


// Handle window resizing
window.addEventListener('resize', () => {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
});
