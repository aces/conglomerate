#include  <bicpl.h>

int  main(
    int   argc,
    char  *argv[] )
{
    Real             min_threshold, max_threshold;
    int              sizes[N_DIMENSIONS], v[MAX_DIMENSIONS], n_dims;
    int              current_label, num_labels, n_neighbours;
    char             *volume_filename, *output_filename;
    Volume           volume, label_volume;
    progress_struct  progress;
    Neighbour_types  connectivity;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( "", &volume_filename ) ||
        !get_string_argument( "", &output_filename ) ||
        !get_real_argument( 0.0, &min_threshold ) ||
        !get_real_argument( 0.0, &max_threshold ) )
    {
        print( "Usage: %s  input.mnc  output.mnc  min_threshold  max_threshold\n",
               argv[0] );
        print( "           [6|26]\n" );

        print( "\n\nCreates a label volume where each connected component has\n" );
        print( "a distinct label number.\n" );
        print( "The connectivity is specified by the last argument as\n" );
        print( "either 6-neighbour or 26-neighbour.\n" );
        return( 1 );
    }

    (void) get_int_argument( 26, &n_neighbours );
    if( n_neighbours == 6 )
        connectivity = FOUR_NEIGHBOURS;
    else if( n_neighbours == 26 )
        connectivity = EIGHT_NEIGHBOURS;
    else
    {
        print( "Connectivity specified must be either 6 or 26.\n" );
        return( 1 );
    }

    if( input_volume( volume_filename, 3, XYZ_dimension_names,
                      NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                      TRUE, &volume, (minc_input_options *) NULL ) != OK )
        return( 1 );

    label_volume = create_label_volume( volume, NC_UNSPECIFIED );

    num_labels = get_volume_voxel_max( label_volume );

    get_volume_sizes( volume, sizes );
    n_dims = get_volume_n_dimensions( volume );

    current_label = 1;

    initialize_progress_report( &progress, FALSE, sizes[n_dims-1],
                                "Filling Regions" );

    BEGIN_ALL_VOXELS( volume, v[0], v[1], v[2], v[3], v[4] )

        if( get_volume_label_data( label_volume, v ) == 0 )
        {
            if( fill_connected_voxels( volume, label_volume, connectivity,
                                       v, 0, 0, current_label,
                                       min_threshold, max_threshold ) &&
                current_label < num_labels-1 )
            {
                ++current_label;
            }
        }

        if( v[n_dims-2] == sizes[n_dims-2] - 1 )
            update_progress_report( &progress, v[n_dims-1] + 1 );

    END_ALL_VOXELS

    terminate_progress_report( &progress );

    print( "Created %d regions.\n", current_label-1 );

    (void) output_volume( output_filename, NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                          label_volume, "cluster_volume",
                          (minc_output_options *) NULL );


    delete_volume( volume );
    delete_volume( label_volume );

    return( 0 );
}