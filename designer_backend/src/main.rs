use axum::{routing::get,
           body::Body,
           http::StatusCode,
           response::{IntoResponse, Response},
           Router,
           Json};
use serde::{Serialize};

#[derive(Serialize)]
struct Cell {
    name: String,
    inferred_type: String,
}


async fn list_available_cells() -> Json<Vec<Cell>>{
    // remember to change this to actually read from the db
    let cells = vec![
        Cell {
            name: "test_b_cell".to_string(),
            inferred_type: "b_cell".to_string(),
        },
        Cell {
            name: "test_t_cell".to_string(),
            inferred_type: "t_cell".to_string()
        }
    ];
    Json(cells)
}

// idea: connect to a database with obs information (including differential expression level)
// grab information for 5 genes and display them as spheres on the surface of the cell
#[tokio::main]
async fn main() {
    let app = Router::new().route("/list_cells", get(list_available_cells));
    println!("[INFO] Running on port 3000");
    // don't forget to change to 0.0.0.0 when we host this
    axum::Server::bind(&"127.0.0.1:3000".parse().expect("[ERROR] Couldn't parse address"))
        .serve(app.into_make_service())
        .await
        .expect("[ERROR] Something went wrong starting the server");
}
