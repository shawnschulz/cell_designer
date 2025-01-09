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
        model.name = 'model';
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

const geometry = new THREE.BoxGeometry( 1, 1, 1 ); 
const material = new THREE.MeshBasicMaterial( {color: 0x00ff00} ); 
const cube = new THREE.Mesh( geometry, material ); 


function placeOnSurface(gltf_object, surface_receptor_gltf_object) {
    // Take a source object, get a random set of vertices, then place
    // the surface receptor gltf_object at the vertex
    scene.add( surface_receptor_gltf_object );
    var verticies = gltf_object.geometry.attributes.position.array
    const vertex = (verticies[0], verticies[1], verticies[2]);
    // test by just placing something at the first vertex

    // Want the world position of the vertex
    gltf_object.localToWorld(vertex);
    surface_receptor_gltf_object.position = vertex;

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
placeOnSurface(model, cube);

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
