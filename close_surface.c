#include  <internal_volume_io.h>
#include  <bicpl.h>

#define  TOLERANCE  1.0e-2

#define  CONNECTIVITY  EIGHT_NEIGHBOURS

#define  BINTREE_FACTOR  0.0

#define  INSIDE_CONVEX_BIT                 128
#define  JUST_INSIDE_CONVEX_HULL_BIT        64
#define  LABEL_MASK         ((JUST_INSIDE_CONVEX_HULL_BIT)-1)

typedef struct
{
    int   voxel[N_DIMENSIONS];
    int   n_voxels;
} convex_boundary_struct;

private  int  remove_just_inside_label(
    Volume   volume,
    int      x,
    int      y,
    int      z );

private  void   label_inside_convex_hull(
    Volume           volume,
    object_struct    *object,
    int              value_to_set );

private  int  label_just_inside_convex_hull(
    Volume                  volume,
    convex_boundary_struct  *boundaries[] );

private  int  get_volume_int_value(
    Volume  volume,
    int     x,
    int     y,
    int     z );

private  BOOLEAN  fill_inside(
    Volume   volume,
    int      x,
    int      y,
    int      z,
    int      label_to_set,
    int      *x_error,
    int      *y_error,
    int      *z_error );

int  main(
    int   argc,
    char  *argv[] )
{
    char                    *input_volume_filename, *input_surface_filename;
    char                    *output_volume_filename;
    int                     i, n_objects, n_convex_boundaries;
    int                     close_threshold;
    int                     x_error, y_error, z_error, label_to_set;
    int                     sizes[N_DIMENSIONS], x, y, z, value;
    convex_boundary_struct  *convex_boundaries;
    STRING                  history;
    File_formats            format;
    Volume                  volume;
    object_struct           **objects;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( "", &input_volume_filename ) ||
        !get_string_argument( "", &input_surface_filename ) ||
        !get_string_argument( "", &output_volume_filename ) )
    {
        print( "Usage: %s  in_volume.mnc  in_surface.obj  out_volume.mnc\n",
               argv[0] );
        print( "   [label_to_set] [close_threshold]\n" );
        return( 1 );
    }

    (void) get_int_argument( 1, &label_to_set );
    (void) get_int_argument( -1, &close_threshold );

    if( input_volume( input_volume_filename, 3, XYZ_dimension_names,
                      NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                      TRUE, &volume, (minc_input_options *) NULL ) != OK )
        return( 1 );

    if( input_graphics_file( input_surface_filename,
                             &format, &n_objects, &objects ) != OK )
        return( 1 );

    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS )
    {
        print( "No objects in %s.\n", input_surface_filename);
        return( 1 );
    }

    label_inside_convex_hull( volume, objects[0], INSIDE_CONVEX_BIT );

    n_convex_boundaries = label_just_inside_convex_hull( volume,
                                                         &convex_boundaries );

    for_less( i, 0, n_convex_boundaries )
    {
        print( "%3d: %d %d %d -- %d voxels\n", i+1,
               convex_boundaries[i].voxel[X],
               convex_boundaries[i].voxel[Y],
               convex_boundaries[i].voxel[Z],
               convex_boundaries[i].n_voxels );
    }

    if( close_threshold > 0 )
    {
        for_less( i, 0, n_convex_boundaries )
        {
            if( convex_boundaries[i].n_voxels >= close_threshold )
            {
                if( remove_just_inside_label( volume,
                        convex_boundaries[i].voxel[X],
                        convex_boundaries[i].voxel[Y],
                        convex_boundaries[i].voxel[Z] ) !=
                        convex_boundaries[i].n_voxels )
                    handle_internal_error( "n voxels" );
            }
        }

        for_less( i, 0, n_convex_boundaries )
        {
            if( convex_boundaries[i].n_voxels >= close_threshold )
            {
                if( fill_inside( volume,
                        convex_boundaries[i].voxel[X],
                        convex_boundaries[i].voxel[Y],
                        convex_boundaries[i].voxel[Z], label_to_set,
                        &x_error, &y_error, &z_error ) )
                {
                    print( "----------------- possible topological hole\n" );
                    print( "%3d: %d %d %d -- leaks through to %d %d %d\n", i+1,
                           convex_boundaries[i].voxel[X],
                           convex_boundaries[i].voxel[Y],
                           convex_boundaries[i].voxel[Z],
                           x_error, y_error, z_error );
                }
            }
        }
    }

    get_volume_sizes( volume, sizes );

    for_less( x, 0, sizes[X] )
    for_less( y, 0, sizes[Y] )
    for_less( z, 0, sizes[Z] )
    {
        value = get_volume_int_value( volume, x, y, z );
        value = value & LABEL_MASK;
        set_volume_real_value( volume, x, y, z, 0, 0, (Real) value );
    }

    (void) strcpy( history, "Inside surface labeled." );

    (void) output_volume( output_volume_filename, NC_UNSPECIFIED,
                          FALSE, 0.0, 0.0, volume, history,
                          (minc_output_options *) NULL );

    return( 0 );
}

private  int  get_volume_int_value(
    Volume  volume,
    int     x,
    int     y,
    int     z )
{
    Real  value;

    value = get_volume_real_value( volume, x, y, z, 0, 0 );

    return( ROUND( value ) );
}

private  void   label_inside_convex_hull(
    Volume           volume,
    object_struct    *object,
    int              value_to_set )
{
    BOOLEAN              inside;
    int                  c, x, y, z, obj_index, n_set;
    int                  sizes[MAX_DIMENSIONS], n_intersects;
    int                  n_points, int_index, next_z;
    Real                 xw, yw, zw, distances[2], limits[2][3];
    Real                 voxel[MAX_DIMENSIONS], max_value, value;
    Real                 boundary_voxel[MAX_DIMENSIONS];
    Point                ray_origin, start_ray, end_ray, *points;
    Point                point_range[2];
    Point                ray_point;
    Vector               ray_direction, offset;
    Real                 **enter_dist, **exit_dist;
    polygons_struct      *polygons;

    polygons = get_polygons_ptr( object );

    max_value = get_volume_real_max( volume );

    if( BINTREE_FACTOR > 0.0 )
    {
        create_polygons_bintree( polygons,
                                 polygons->n_items * BINTREE_FACTOR + 1);
    }

    n_points = polygons->n_points;
    points = polygons->points;

    get_range_points( n_points, points, &point_range[0], &point_range[1] );

    for_less( x, 0, 2 )
    for_less( y, 0, 2 )
    for_less( z, 0, 2 )
    {
        convert_world_to_voxel( volume, Point_x(point_range[x]),
                                Point_y(point_range[y]),
                                Point_z(point_range[z]), voxel );
        for_less( c, 0, N_DIMENSIONS )
        {
            if( x == 0 && y == 0 && z == 0 || voxel[c] < limits[0][c] )
            {
                limits[0][c] = voxel[c];
            }
            if( x == 0 && y == 0 && z == 0 || voxel[c] > limits[1][c] )
            {
                limits[1][c] = voxel[c];
            }
        }

    }

    for_less( c, 0, N_DIMENSIONS )
    {
        limits[0][c] -= 100.0;
        limits[1][c] += 100.0;
    }

    get_volume_sizes( volume, sizes );

    n_set = 0;

    ALLOC2D( enter_dist, sizes[X], sizes[Y] );
    ALLOC2D( exit_dist, sizes[X], sizes[Y] );

    for_less( x, 0, sizes[X] )
    {
        voxel[X] = (Real) x;
        for_less( y, 0, sizes[Y] )
        {
            voxel[Y] = (Real) y;

            voxel[Z] = limits[0][Z];
            convert_voxel_to_world( volume, voxel, &xw, &yw, &zw );
            fill_Point( start_ray, xw, yw, zw );

            voxel[Z] = limits[1][Z];
            convert_voxel_to_world( volume, voxel, &xw, &yw, &zw );
            fill_Point( end_ray, xw, yw, zw );

            ray_origin = end_ray;
            SUB_VECTORS( ray_direction, start_ray, ray_origin );
            NORMALIZE_VECTOR( ray_direction, ray_direction );

            enter_dist[x][y] = -1.0;
            exit_dist[x][y] = -1.0;

            if( intersect_ray_with_object( &ray_origin, &ray_direction,
                                           object, &obj_index,
                                           &distances[1], NULL ) == 0 )
                continue;
   
            SUB_VECTORS( offset, end_ray, start_ray );
            distances[1] = MAGNITUDE( offset ) - distances[1];
            exit_dist[x][y] = distances[1];

            ray_origin = start_ray;
            SUB_VECTORS( ray_direction, end_ray, ray_origin );
            NORMALIZE_VECTOR( ray_direction, ray_direction );

            if( intersect_ray_with_object( &ray_origin, &ray_direction,
                                           object, &obj_index,
                                           &distances[0], NULL ) == 0 )
            {
                print( "ray distance error\n" );
            }

            enter_dist[x][y] = distances[0];

            n_intersects = 2;

            inside = FALSE;
            int_index = -1;
            next_z = -1;

            for_less( z, 0, sizes[Z] )
            {
                while( next_z <= z )
                {
                    ++int_index;
                    if( int_index < n_intersects )
                    {
                        GET_POINT_ON_RAY( ray_point, ray_origin, ray_direction,
                                          distances[int_index] );
                        convert_world_to_voxel( volume, 
                                                Point_x(ray_point),
                                                Point_y(ray_point),
                                                Point_z(ray_point),
                                                boundary_voxel );
                        next_z = CEILING( boundary_voxel[Z] );
                        inside = ((int_index % 2) == 1);
                    }
                    else
                    {
                        next_z = sizes[Z];
                        inside = FALSE;
                    }
                }

                if( inside )
                {
                    value = get_volume_real_value( volume, x, y, z, 0, 0);
                    value = (int) value | value_to_set;
                    if( value > max_value )
                        value = max_value;
                    set_volume_real_value( volume, x, y, z, 0, 0, value);
                    ++n_set;
                }
            }
        }
    }

    print( "Set %d out of %d\n", n_set, sizes[X] * sizes[Y] * sizes[Z] );

    for_less( x, 1, sizes[X]-1 )
    {
        int      dx, dy, dir;
        BOOLEAN  error;

        for_less( y, 1, sizes[Y]-1 )
        {
            error = FALSE;

            for_less( dir, 0, 4 )
            {
                switch( dir )
                {
                case 0:  dx = 1;  dy = 0;  break;
                case 1:  dx = 1;  dy = 1;  break;
                case 2:  dx = 0;  dy = 1;  break;
                case 3:  dx = -1;  dy = 1;  break;
                }

                if( enter_dist[x-dx][y-dy] >= 0.0 && enter_dist[x+dx][y+dy] >= 0.0&&
                    (enter_dist[x][y] < 0.0 || enter_dist[x][y] - TOLERANCE >
                     (enter_dist[x-dx][y-dy] + enter_dist[x+dx][y+dy])/2.0) )
                    error = TRUE;

                if( exit_dist[x-dx][y-dy] >= 0.0 && exit_dist[x+dx][y+dy] >= 0.0&&
                    (exit_dist[x][y] < 0.0 || exit_dist[x][y] + TOLERANCE <
                     (exit_dist[x-dx][y-dy] + exit_dist[x+dx][y+dy])/2.0) )
                    error = TRUE;
            }

            if( error )
            {
                print( "%d %d: ", x, y );
                handle_internal_error( "enter_dist" );
            }
        }
    }

    FREE2D( enter_dist );
    FREE2D( exit_dist );
}

private  BOOLEAN  is_on_convex_boundary(
    Volume   volume,
    int      sizes[],
    int      x,
    int      y,
    int      z )
{
    int   i, n_dirs, *dx, *dy, *dz, tx, ty, tz, value;

    value = get_volume_int_value( volume, x, y, z );
    if( (value & LABEL_MASK) != 0 ||
        ((value & INSIDE_CONVEX_BIT) == 0 ||
         (value & JUST_INSIDE_CONVEX_HULL_BIT) != 0) )
        return( FALSE );

    n_dirs = get_3D_neighbour_directions( CONNECTIVITY, &dx, &dy, &dz );

    for_less( i, 0, n_dirs )
    {
        tx = x + dx[i];
        ty = y + dy[i];
        tz = z + dz[i];

        if( tx < 0 || tx >= sizes[X] ||
            ty < 0 || ty >= sizes[Y] ||
            tz < 0 || tz >= sizes[Z] ||
            (get_volume_int_value( volume, tx, ty, tz ) & INSIDE_CONVEX_BIT)
                                                            == 0 )
        {
            return( TRUE );
        }
    }

    return( FALSE );
}

typedef struct
{
    short  x, y, z;
} xyz_struct;

private  int  expand_convex_boundary(
    Volume   volume,
    int      sizes[],
    int      x,
    int      y,
    int      z )
{
    int                          i, n_dirs, *dx, *dy, *dz, tx, ty, tz, value;
    int                          n_voxels;
    xyz_struct                   xyz;
    QUEUE_STRUCT( xyz_struct )   queue;

    INITIALIZE_QUEUE( queue );

    value = get_volume_int_value( volume, x, y, z );
    value |= JUST_INSIDE_CONVEX_HULL_BIT;
    set_volume_real_value( volume, x, y, z, 0, 0, (Real) value );
    xyz.x = (short) x;
    xyz.y = (short) y;
    xyz.z = (short) z;

    n_dirs = get_3D_neighbour_directions( CONNECTIVITY, &dx, &dy, &dz );

    INSERT_IN_QUEUE( queue, xyz );

    n_voxels = 1;

    while( !IS_QUEUE_EMPTY( queue ) )
    {
        REMOVE_FROM_QUEUE( queue, xyz );

        x = (int) xyz.x;
        y = (int) xyz.y;
        z = (int) xyz.z;

        for_less( i, 0, n_dirs )
        {
            tx = x + dx[i];
            ty = y + dy[i];
            tz = z + dz[i];

            if( tx >= 0 && tx < sizes[X] &&
                ty >= 0 && ty < sizes[Y] &&
                tz >= 0 && tz < sizes[Z] &&
                is_on_convex_boundary( volume, sizes, tx, ty, tz ) )
            {
                value = get_volume_int_value( volume, tx, ty, tz );
                value |= JUST_INSIDE_CONVEX_HULL_BIT;
                set_volume_real_value( volume, tx, ty, tz, 0, 0, (Real) value );
                xyz.x = (short) tx;
                xyz.y = (short) ty;
                xyz.z = (short) tz;

                INSERT_IN_QUEUE( queue, xyz );
                ++n_voxels;
            }
        }
    }

    DELETE_QUEUE( queue );

    return( n_voxels );
}

private  int  label_just_inside_convex_hull(
    Volume                  volume,
    convex_boundary_struct  *boundaries[] )
{
    int                     sizes[N_DIMENSIONS], x, y, z, n_boundaries;
    convex_boundary_struct  bound;

    get_volume_sizes( volume, sizes );
    n_boundaries = 0;

    for_less( x, 0, sizes[X] )
    for_less( y, 0, sizes[Y] )
    for_less( z, 0, sizes[Z] )
    {
        if( is_on_convex_boundary( volume, sizes, x, y, z ) )
        {
            bound.voxel[X] = x;
            bound.voxel[Y] = y;
            bound.voxel[Z] = z;
            bound.n_voxels = expand_convex_boundary( volume, sizes, x, y, z );
            ADD_ELEMENT_TO_ARRAY( *boundaries, n_boundaries, bound,
                                  DEFAULT_CHUNK_SIZE );
        }
    }

    return( n_boundaries );
}

private  int  remove_just_inside_label(
    Volume   volume,
    int      x,
    int      y,
    int      z )
{
    int                          i, n_dirs, *dx, *dy, *dz, tx, ty, tz, value;
    int                          n_voxels, sizes[N_DIMENSIONS];
    xyz_struct                   xyz;
    QUEUE_STRUCT( xyz_struct )   queue;

    get_volume_sizes( volume, sizes );

    INITIALIZE_QUEUE( queue );

    value = get_volume_int_value( volume, x, y, z );
    value -= JUST_INSIDE_CONVEX_HULL_BIT;
    set_volume_real_value( volume, x, y, z, 0, 0, (Real) value );
    xyz.x = (short) x;
    xyz.y = (short) y;
    xyz.z = (short) z;

    n_dirs = get_3D_neighbour_directions( CONNECTIVITY, &dx, &dy, &dz );

    INSERT_IN_QUEUE( queue, xyz );

    n_voxels = 1;

    while( !IS_QUEUE_EMPTY( queue ) )
    {
        REMOVE_FROM_QUEUE( queue, xyz );

        x = (int) xyz.x;
        y = (int) xyz.y;
        z = (int) xyz.z;

        for_less( i, 0, n_dirs )
        {
            tx = x + dx[i];
            ty = y + dy[i];
            tz = z + dz[i];

            if( tx >= 0 && tx < sizes[X] &&
                ty >= 0 && ty < sizes[Y] &&
                tz >= 0 && tz < sizes[Z] &&
                (get_volume_int_value( volume, tx, ty, tz ) &
                    JUST_INSIDE_CONVEX_HULL_BIT) != 0 )
            {
                value = get_volume_int_value( volume, tx, ty, tz );
                value -= JUST_INSIDE_CONVEX_HULL_BIT;
                set_volume_real_value( volume, tx, ty, tz, 0, 0, (Real) value );
                xyz.x = (short) tx;
                xyz.y = (short) ty;
                xyz.z = (short) tz;

                INSERT_IN_QUEUE( queue, xyz );
                ++n_voxels;
            }
        }
    }

    DELETE_QUEUE( queue );

    return( n_voxels );
}

private  BOOLEAN  fill_inside(
    Volume   volume,
    int      x,
    int      y,
    int      z,
    int      label_to_set,
    int      *x_error,
    int      *y_error,
    int      *z_error )
{
    BOOLEAN                      error;
    int                          i, n_dirs, *dx, *dy, *dz, tx, ty, tz, value;
    int                          sizes[N_DIMENSIONS];
    xyz_struct                   xyz;
    QUEUE_STRUCT( xyz_struct )   queue;

    error = FALSE;

    get_volume_sizes( volume, sizes );

    INITIALIZE_QUEUE( queue );

    value = get_volume_int_value( volume, x, y, z );
    value |= label_to_set;
    set_volume_real_value( volume, x, y, z, 0, 0, (Real) value );
    xyz.x = (short) x;
    xyz.y = (short) y;
    xyz.z = (short) z;

    n_dirs = get_3D_neighbour_directions( CONNECTIVITY, &dx, &dy, &dz );

    INSERT_IN_QUEUE( queue, xyz );

    while( !IS_QUEUE_EMPTY( queue ) )
    {
        REMOVE_FROM_QUEUE( queue, xyz );

        x = (int) xyz.x;
        y = (int) xyz.y;
        z = (int) xyz.z;

        for_less( i, 0, n_dirs )
        {
            tx = x + dx[i];
            ty = y + dy[i];
            tz = z + dz[i];

            if( tx >= 0 && tx < sizes[X] &&
                ty >= 0 && ty < sizes[Y] &&
                tz >= 0 && tz < sizes[Z] )
            {
                value = get_volume_int_value( volume, tx, ty, tz );

                if( (value & LABEL_MASK) == 0 &&
                    (value & INSIDE_CONVEX_BIT) != 0 )
                {
                    if( (value & JUST_INSIDE_CONVEX_HULL_BIT) != 0 &&
                        !error )
                    {
                        error = TRUE;
                        *x_error = tx;
                        *y_error = ty;
                        *z_error = tz;
                    }

                    value |= label_to_set;
                    set_volume_real_value( volume, tx, ty, tz, 0, 0,
                                           (Real) value );
                    xyz.x = (short) tx;
                    xyz.y = (short) ty;
                    xyz.z = (short) tz;
    
                    INSERT_IN_QUEUE( queue, xyz );
                }
            }
        }
    }

    DELETE_QUEUE( queue );

    return( error );
}

