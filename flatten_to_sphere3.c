#include  <internal_volume_io.h>
#include  <bicpl.h>

private  void  flatten_polygons(
    polygons_struct  *polygons,
    Point            init_points[],
    Real             sphere_weight,
    int              n_iters );

int  main(
    int    argc,
    char   *argv[] )
{
    STRING               src_filename, dest_filename, initial_filename;
    int                  n_objects, n_i_objects, n_iters;
    File_formats         format;
    object_struct        **object_list, **i_object_list;
    polygons_struct      *polygons, *init_polygons;
    Point                *init_points;
    Real                 sphere_weight;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &src_filename ) ||
        !get_string_argument( NULL, &dest_filename ) )
    {
        print_error( "Usage: %s  input.obj output.obj [n_iters]\n",
                     argv[0] );
        return( 1 );
    }

    (void) get_real_argument( 1.0, &sphere_weight );
    (void) get_int_argument( 100, &n_iters );

    if( input_graphics_file( src_filename, &format, &n_objects,
                             &object_list ) != OK || n_objects != 1 ||
        get_object_type(object_list[0]) != POLYGONS )
        return( 1 );

    polygons = get_polygons_ptr( object_list[0] );

    if( get_string_argument( NULL, &initial_filename ) )
    {
        if( input_graphics_file( initial_filename, &format, &n_i_objects,
                                 &i_object_list ) != OK || n_i_objects != 1 ||
            get_object_type(i_object_list[0]) != POLYGONS )
            return( 1 );

        init_polygons = get_polygons_ptr( i_object_list[0] );
        init_points = init_polygons->points;
        ALLOC( init_polygons->points, 1 );
        delete_object_list( n_i_objects, i_object_list );
    }
    else
    {
        init_points = polygons->points;
    }

    flatten_polygons( polygons, init_points, sphere_weight, n_iters );

    if( output_graphics_file( dest_filename, format, 1, object_list ) != OK )
        print_error( "Error outputting: %s\n", dest_filename );

    return( 0 );
}

private  Real  evaluate_fit(
    int     n_parameters,
    Real    parameters[],
    Real    distances[],
    int     n_neighbours[],
    int     *neighbours[],
    Real    sphere_weight )
{
    int    p, n_points, ind, n, neigh;
    Real   fit, xc, yc, zc, x, y, z, dx, dy, dz, dist, diff, act_dist;

    fit = 0.0;
    n_points = n_parameters / 3;

    ind = 0;
    for_less( p, 0, n_points )
    {
        for_less( n, 0, n_neighbours[p] )
        {
            neigh = neighbours[p][n];
            if( neigh < p )
                continue;

            dist = distances[ind];
            ++ind;

            dx = parameters[IJ(p,0,3)] - parameters[IJ(neigh,0,3)];
            dy = parameters[IJ(p,1,3)] - parameters[IJ(neigh,1,3)];
            dz = parameters[IJ(p,2,3)] - parameters[IJ(neigh,2,3)];
            act_dist = dx * dx + dy * dy + dz * dz;
            diff = dist - act_dist;
            fit += diff * diff;
        }
    }

    for_less( p, 0, n_points )
    {
        xc = 0.0;
        yc = 0.0;
        zc = 0.0;
        for_less( n, 0, n_neighbours[p] )
        {
            neigh = neighbours[p][n];
            xc += parameters[IJ(neigh,0,3)];
            yc += parameters[IJ(neigh,1,3)];
            zc += parameters[IJ(neigh,2,3)];
        }

        xc /= (Real) n_neighbours[p];
        yc /= (Real) n_neighbours[p];
        zc /= (Real) n_neighbours[p];

        x = parameters[IJ(p,0,3)];
        y = parameters[IJ(p,1,3)];
        z = parameters[IJ(p,2,3)];

        dx = x - xc;
        dy = y - yc;
        dz = z - zc;

        fit += sphere_weight * (dx * dx + dy * dy + dz * dz);
    }

    return( fit );
}

private  void  evaluate_fit_derivative(
    int     n_parameters,
    Real    parameters[],
    Real    distances[],
    int     n_neighbours[],
    int     *neighbours[],
    Real    sphere_weight,
    Real    deriv[] )
{
    int    p, n_points, ind, n, neigh, dim;
    Real   dx, dy, dz, dist, diff, act_dist, factor;
    Real   x1, y1, z1, x2, y2, z2, xc;

    for_less( p, 0, n_parameters )
        deriv[p] = 0.0;

    n_points = n_parameters / 3;

    ind = 0;
    for_less( p, 0, n_points )
    {
        for_less( n, 0, n_neighbours[p] )
        {
            neigh = neighbours[p][n];
            if( neigh < p )
                continue;

            dist = distances[ind];
            ++ind;

            x1 = parameters[IJ(p,0,3)];
            y1 = parameters[IJ(p,1,3)];
            z1 = parameters[IJ(p,2,3)];
            x2 = parameters[IJ(neigh,0,3)];
            y2 = parameters[IJ(neigh,1,3)];
            z2 = parameters[IJ(neigh,2,3)];
            dx = x1 - x2;
            dy = y1 - y2;
            dz = z1 - z2;
            act_dist = dx * dx + dy * dy + dz * dz;
            factor = act_dist - dist;
            deriv[IJ(p,0,3)] += 2.0 * (x1 - x2) * factor;
            deriv[IJ(p,1,3)] += 2.0 * (y1 - y2) * factor;
            deriv[IJ(p,2,3)] += 2.0 * (z1 - z2) * factor;
            deriv[IJ(neigh,0,3)] += 2.0 * (x2 - x1) * factor;
            deriv[IJ(neigh,1,3)] += 2.0 * (y2 - y1) * factor;
            deriv[IJ(neigh,2,3)] += 2.0 * (z2 - z1) * factor;
        }
    }

    for_less( p, 0, n_points )
    {
        for_less( dim, 0, N_DIMENSIONS )
        {
            xc = 0.0;
            for_less( n, 0, n_neighbours[p] )
            {
                neigh = neighbours[p][n];
                xc += parameters[IJ(neigh,dim,3)];
            }

            xc /= (Real) n_neighbours[p];
            diff = xc - parameters[IJ(p,dim,3)];
            deriv[IJ(p,dim,3)] += sphere_weight * -2.0 * diff * 1.0;
            for_less( n, 0, n_neighbours[p] )
            {
                neigh = neighbours[p][n];
                deriv[IJ(neigh,dim,3)] += sphere_weight * 2.0 * diff /
                                          (Real) n_neighbours[p];
            }
        }
    }
}

private  void  evaluate_fit_along_line(
    int     n_parameters,
    Real    parameters[],
    Real    delta[],
    Real    distances[],
    int     n_neighbours[],
    int     *neighbours[],
    Real    sphere_weight,
    Real    coefs[] )
{
    int    p, n_points, ind, n, neigh, dim;
    Real   dx, dy, dz, dist;
    Real   x, y, z, weight;
    Real   line_coefs[3];

    for_less( p, 0, 5 )
        coefs[p] = 0.0;

    n_points = n_parameters / 3;

    ind = 0;
    for_less( p, 0, n_points )
    {
        for_less( n, 0, n_neighbours[p] )
        {
            neigh = neighbours[p][n];
            if( neigh < p )
                continue;

            dist = distances[ind];
            ++ind;

            x = parameters[IJ(p,0,3)] - parameters[IJ(neigh,0,3)];
            y = parameters[IJ(p,1,3)] - parameters[IJ(neigh,1,3)];
            z = parameters[IJ(p,2,3)] - parameters[IJ(neigh,2,3)];
            dx = delta[IJ(p,0,3)] - delta[IJ(neigh,0,3)];
            dy = delta[IJ(p,1,3)] - delta[IJ(neigh,1,3)];
            dz = delta[IJ(p,2,3)] - delta[IJ(neigh,2,3)];

            line_coefs[0] = x * x + y * y + z * z - dist;
            line_coefs[1] = 2.0 * (x * dx + y * dy + z * dz);
            line_coefs[2] = dx * dx + dy * dy + dz * dz;

            coefs[0] += line_coefs[0] * line_coefs[0];
            coefs[1] += 2.0 * line_coefs[1] * line_coefs[0];
            coefs[2] += 2.0 * (line_coefs[2] * line_coefs[0] +
                               line_coefs[1] * line_coefs[1]);
            coefs[3] += 2.0 * line_coefs[2] * line_coefs[1];
            coefs[4] += line_coefs[2] * line_coefs[2];
        }
    }

    for_less( p, 0, n_points )
    {
        for_less( dim, 0, N_DIMENSIONS )
        {
            line_coefs[0] = -parameters[IJ(p,dim,3)];
            line_coefs[1] = -delta[IJ(p,dim,3)];
            weight = 1.0 / (Real) n_neighbours[p];

            for_less( n, 0, n_neighbours[p] )
            {
                neigh = neighbours[p][n];
                line_coefs[0] += parameters[IJ(neigh,dim,3)] * weight;
                line_coefs[1] += delta[IJ(neigh,dim,3)] * weight;
            }

            coefs[0] += sphere_weight * line_coefs[0] * line_coefs[0];
            coefs[1] += sphere_weight * 2.0 * line_coefs[0] * line_coefs[1];
            coefs[2] += sphere_weight * line_coefs[1] * line_coefs[1];
        }
    }
}

private  void  minimize_along_line(
    int     n_parameters,
    Real    parameters[],
    Real    delta[],
    Real    distances[],
    int     n_neighbours[],
    int     *neighbours[],
    Real    sphere_weight )
{
    int    p, s, n_solutions, best_index;
    Real   coefs[5], deriv[4], t, fit, best_fit, solutions[3];
#ifdef   DEBUG
#define  DEBUG
    Real   test_fit, *test;

    ALLOC( test, n_parameters );
#endif

    evaluate_fit_along_line( n_parameters, parameters, delta, distances,
                             n_neighbours, neighbours, sphere_weight,
                             coefs );

    for_less( p, 0, 4 )
        deriv[p] = (Real) (p+1) * coefs[p+1];

    n_solutions = solve_cubic( deriv[3], deriv[2], deriv[1], deriv[0],
                               solutions );

    if( n_solutions == 0 )
        print( "N solutions = 0\n" );

    best_fit = coefs[0];
    best_index = -1;

#ifdef DEBUG
    test_fit = evaluate_fit( n_parameters, parameters, distances,
                             n_neighbours, neighbours, sphere_weight );
    print( "## %g %g\n", coefs[0], test_fit );
#endif

    for_less( s, 0, n_solutions )
    {
        t = solutions[s];

        fit = coefs[0] + t * (coefs[1] + t * (coefs[2] + t * (coefs[3] +
              t * coefs[4])));

#ifdef DEBUG
        for_less( p, 0, n_parameters )
            test[p] = parameters[p] + t * delta[p];

        test_fit = evaluate_fit( n_parameters, test, distances,
                                 n_neighbours, neighbours, sphere_weight );

        print( "%g %g\n", fit, test_fit );
#endif

        if( fit < best_fit )
        {
            best_fit = fit;
            best_index = s;
        }
    }

    if( best_index >= 0 )
    {
        t = solutions[best_index];

        for_less( p, 0, n_parameters )
            parameters[p] = parameters[p] + t * delta[p];
    }

#ifdef DEBUG
    FREE( test );
#endif
}

private  void  flatten_polygons(
    polygons_struct  *polygons,
    Point            init_points[],
    Real             sphere_weight,
    int              n_iters )
{
    int              p, n, point, *n_neighbours, **neighbours;
    int              n_parameters, total_neighbours;
    Real             gg, dgg, gam, current_time, last_update_time;
    Real             *parameters, *g, *h, *xi, fit, *unit_dir;
    Real             *distances, len;
    int              iter, ind, update_rate;

    create_polygon_point_neighbours( polygons, FALSE, &n_neighbours,
                                     &neighbours, NULL, NULL );

    total_neighbours = 0;
    for_less( point, 0, polygons->n_points )
        total_neighbours += n_neighbours[point];
    total_neighbours /= 2;

    ALLOC( distances, total_neighbours );
    ind = 0;

    for_less( point, 0, polygons->n_points )
    {
        for_less( n, 0, n_neighbours[point] )
        {
            if( neighbours[point][n] < point )
                continue;
            distances[ind] = sq_distance_between_points(
                                   &polygons->points[point],
                                   &polygons->points[neighbours[point][n]] );
            ++ind;
        }
    }

    n_parameters = 3 * polygons->n_points;
    ALLOC( parameters, n_parameters );
    ALLOC( g, n_parameters );
    ALLOC( h, n_parameters );
    ALLOC( xi, n_parameters );
    ALLOC( unit_dir, n_parameters );

    for_less( point, 0, polygons->n_points )
    {
        parameters[IJ(point,0,3)] = RPoint_x(init_points[point] );
        parameters[IJ(point,1,3)] = RPoint_y(init_points[point] );
        parameters[IJ(point,2,3)] = RPoint_z(init_points[point] );
    }

    sphere_weight *= (Real) total_neighbours / (Real) polygons->n_points;

    fit = evaluate_fit( n_parameters, parameters, distances,
                        n_neighbours, neighbours, sphere_weight );

    print( "Initial  %g\n", fit );
    (void) flush_file( stdout );

    evaluate_fit_derivative( n_parameters, parameters, distances,
                             n_neighbours, neighbours, sphere_weight, xi );

    for_less( p, 0, n_parameters )
    {
        g[p] = -xi[p];
        h[p] = g[p];
        xi[p] = g[p];
    }

    update_rate = 1;
    last_update_time = current_cpu_seconds();

    for_less( iter, 0, n_iters )
    {
        len = 0.0;
        for_less( p, 0, n_parameters )
            len += xi[p] * xi[p];

        len = sqrt( len );
        for_less( p, 0, n_parameters )
            unit_dir[p] = xi[p] / len;

        minimize_along_line( n_parameters, parameters, unit_dir, distances,
                             n_neighbours, neighbours, sphere_weight );

        if( ((iter+1) % update_rate) == 0 || iter == n_iters - 1 )
        {
            fit = evaluate_fit( n_parameters, parameters, distances,
                                n_neighbours, neighbours, sphere_weight );

            print( "%d: %g\n", iter+1, fit );
            (void) flush_file( stdout );
            current_time = current_cpu_seconds();
            if( current_time - last_update_time < 1.0 )
                update_rate *= 10;
            last_update_time = current_time;
        }

        evaluate_fit_derivative( n_parameters, parameters, distances,
                                 n_neighbours, neighbours, sphere_weight, xi );

        gg = 0.0;
        dgg = 0.0;
        for_less( p, 0, n_parameters )
        {
            gg += g[p] * g[p];
            dgg += (xi[p] + g[p]) * xi[p];
/*
            dgg += xi[p] * xi[p];
*/
        }

        if( len == 0.0 )
            break;

        gam = dgg / gg;

        for_less( p, 0, n_parameters )
        {
            g[p] = -xi[p];
            h[p] = g[p] + gam * h[p];
            xi[p] = h[p];
        }
    }

    delete_polygon_point_neighbours( polygons, n_neighbours, neighbours,
                                     NULL, NULL );

    for_less( point, 0, polygons->n_points )
    {
        fill_Point( polygons->points[point],
                    parameters[IJ(point,0,3)],
                    parameters[IJ(point,1,3)],
                    parameters[IJ(point,2,3)] );
    }

    FREE( distances );
    FREE( parameters );
    FREE( xi );
    FREE( g );
    FREE( h );
    FREE( unit_dir );
}
