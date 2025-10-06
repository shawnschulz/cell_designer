use axum::{routing::get, Router};

// idea: connect to a database with obs information (including differential expression level)
// grab information for 5 genes and display them as spheres on the surface of the cell
#[tokio::main]
async fn main() {
    let app = Router::new().route("/", get(|| async {"test_cell"}));
    println!("[INFO] Running on port 3000");
    // don't forget to change to 0.0.0.0 when we host this
    axum::Server::bind(&"localhost:3000".parse().unwrap())
        .serve(app.into_make_service())
        .await()
        .expect("[ERROR] Server failed to start")

}
