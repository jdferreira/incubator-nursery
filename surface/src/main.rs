extern crate image;

mod surface;

const N: u32 = 300;

fn main() {
    
    let mut surface = surface::Surface::new();
    
    surface.add_particle(surface::SurfacePoint::new_degree(-45.0, -45.0), -0.5, 1.0);
    surface.add_particle(surface::SurfacePoint::new_degree(150.0, -80.0), -0.2, 0.1);
    surface.add_particle(surface::SurfacePoint::new_degree( 79.0,  80.0),  0.7, 0.5);
    surface.add_particle(surface::SurfacePoint::new_degree( 45.0,   0.0), -0.7, 0.8);
    
    let p1 = surface.add_point(surface::SurfacePoint::new_degree( 10.0, -80.0));
    let p2 = surface.add_point(surface::SurfacePoint::new_degree( 10.0, -10.0));
    let p3 = surface.add_point(surface::SurfacePoint::new_degree( 90.0,   0.0));
    let p4 = surface.add_point(surface::SurfacePoint::new_degree(180.0,  10.0));
    let p5 = surface.add_point(surface::SurfacePoint::new_degree(270.0,   0.0));
    let p6 = surface.add_point(surface::SurfacePoint::new_degree(350.0,  80.0));
        
    surface.add_edge(p1, p2);
    surface.add_edge(p1, p3);
    surface.add_edge(p1, p4);
    surface.add_edge(p1, p5);
    
    surface.add_edge(p6, p2);
    surface.add_edge(p6, p3);
    surface.add_edge(p6, p4);
    surface.add_edge(p6, p5);
    
    surface.add_edge(p2, p3);
    surface.add_edge(p3, p4);
    surface.add_edge(p4, p5);
    surface.add_edge(p5, p2);
    
    let background = surface::get_background(&surface, 600);
    surface::draw(&surface, &background, "img/map0.png");
    
    for i in 1..(N + 1) {
        surface.step(0.1);
        if i % 10 == 0 {
            println!("{:?}", i);
            surface::draw(&surface, &background, &format!("img/map{}.png", i));
        }
    }
    
    surface::draw(&surface, &background, &format!("img/map{}.png", N));
}
