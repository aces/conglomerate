/*
*---------------------------------------------------------------------------
*@COPYRIGHT :
*             Copyright 2005, Jonathan Harlap
*             McConnell Brain Imaging Centre,
*             Montreal Neurological Institute, McGill University.
*             Permission to use, copy, modify, and distribute this
*             software and its documentation for any purpose and without
*             fee is hereby granted, provided that the above copyright
*             notice appear in all copies.  The author and McGill University
*             make no representations about the suitability of this
*             software for any purpose.  It is provided "as is" without
*             express or implied warranty.
*---------------------------------------------------------------------------- 
*$Revision: 1.5 $
*$Author: clepage $
*$Date: 2018/11/29 17:58:46 $
*---------------------------------------------------------------------------
*
* print_world_values <minclist> <coordlist> <outputfile>
* 
* Reads the list of minc files in <minclist> (note that to handle
* glim_image files, it only reads the first word per line).  Reads the
* list of world coordinates in <coordlist> assuming that each line
* contains space separated x, y, and z coordinates -
* implicitly this means that this program only deals with 3D volumes.
* Finally, for each minc it reads the real value of each coordinate
* specified in the coordlist and writes it out as a tab separated list
* into the output file.
*/

#include  <stdio.h>
#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>

void usage(char* progname) {
        printf("Usage: %s <glimfile> <coordlist file> <outputfile>\n\nNote that the glimfile is the list of minc files, using the first column of the file as the list of minc files.\nThe coordlist file contains one world coordinate in the form 'x y z' per line.\n", progname);
}

#define EVAL_NEAREST -1
#define EVAL_LINEAR 0
#define EVAL_CUBIC 2

int interpolant = EVAL_LINEAR;

int  main(
                         int   argc,
                         char  *argv[] )
{
        STRING     glim_filename, coordlist_filename, output_filename;
        char       cur_minc[255];
        double     *x, *y, *z;
        double     value;
        Volume     volume;
        double     curx, cury, curz;
        int        voxx, voxy, voxz, sizes[MAX_DIMENSIONS];
        Real       voxel[MAX_DIMENSIONS];
        int        i, r, keep_looping, n_coords;
        FILE*      coordfile;
        FILE*      glimfile;
        FILE*      outputfile;
        
        initialize_argument_processing( argc, argv );

        if( !get_string_argument( "", &glim_filename ) ) {
          usage(argv[0]);
          return( 1 );
        }
 
        if( !get_string_argument( "", &coordlist_filename ) ) {
          usage(argv[0]);
          return( 1 );
        }
 
        if( !get_string_argument( "", &output_filename ) ) {
          usage(argv[0]);
          return( 1 );
        }
 
        /* first pass: count the number of coordinates in the list */
        coordfile = fopen(coordlist_filename, "r+t");
        n_coords = 0;
        keep_looping = 1;
        while(keep_looping) {
          if(fscanf(coordfile, "%lf%lf%lf", &curx, &cury, &curz) != 3) {
            keep_looping = 0;
          } else {
            n_coords++;
          }
        }
        fclose(coordfile);

        x = (double*)malloc( n_coords * sizeof( double ) );
        y = (double*)malloc( n_coords * sizeof( double ) );
        z = (double*)malloc( n_coords * sizeof( double ) );
        if( !x || !y || !z ) {
          printf( "Failed to allocate memory for %d coordinates.\n", 
                  n_coords );
          return(1);
        }

        /* read coordlist into x/y/z arrays */
        coordfile = fopen(coordlist_filename, "r+t");
        i = 0;
        keep_looping = 1;
        while(keep_looping && i < n_coords) {
          curx = 0;
          cury = 0;
          curz = 0;
          if(fscanf(coordfile, "%lf%lf%lf", &curx, &cury, &curz) != 3) {
            keep_looping = 0;
          } else {
            x[i] = curx;
            y[i] = cury;
            z[i] = curz;
            ++i;
          }
        }
        fclose(coordfile);

        printf("There are %d coords.\n", n_coords);
        for(i = 0; i < n_coords; ++i) {
                printf("%d: %lf %lf %lf\n", i, x[i], y[i], z[i]);
        }
 
        /* read the glim file to get the minc filenames */
        keep_looping = 1;
        glimfile = fopen(glim_filename, "r+t");
        if( glimfile == NULL ) {
          printf( "Error: Cannot open glim file %s\n", glim_filename );
        }
        outputfile = fopen(output_filename, "w+t");
        if( outputfile == NULL ) {
          printf( "Error: Cannot open output file %s\n", output_filename );
        }

        /* print out the x,y,z of each coordinate */
        fprintf(outputfile, "x\t");
        for(i = 0; i < n_coords; ++i)
          fprintf(outputfile, "%lf\t", x[i]);
        fprintf(outputfile, "\ny\t");
        for(i = 0; i < n_coords; ++i)
          fprintf(outputfile, "%lf\t", y[i]);
        fprintf(outputfile, "\nz\t");
        for(i = 0; i < n_coords; ++i)
          fprintf(outputfile, "%lf\t", z[i]);
        fprintf(outputfile, "\n");
          
        while(keep_looping) {
                if(fscanf(glimfile, " %[^ \t\n] %*1[\n]", &cur_minc) > 0) {

                        if( input_volume( cur_minc, 3, XYZ_dimension_names,
                                                                        NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                                                                        TRUE, &volume, (minc_input_options *) NULL ) != OK ) {
                                printf("Failed to read %s\n", cur_minc);
                                return( 1 );
                        }
                        
                        fprintf(outputfile, "%s\t", cur_minc);

                        for(i = 0; i < n_coords ; ++i) {

                                get_volume_sizes( volume, sizes );

                                convert_world_to_voxel(volume, x[i], y[i], z[i], voxel);

                                voxx = FLOOR( voxel[0] );
                                voxy = FLOOR( voxel[1] );
                                voxz = FLOOR( voxel[2] );

                                printf("Reading %lf %lf %lf (%d %d %d) from %s\n", x[i], y[i], z[i], voxx, voxy, voxz, cur_minc);

                                if( voxx >= 0 && voxx < sizes[0] &&
                                         voxy >= 0 && voxy < sizes[1] &&
                                         voxz >= 0 && voxz < sizes[2] )
                                        {

                                                evaluate_volume_in_world( volume, x[i], y[i], z[i], interpolant, FALSE, 0.0, &value, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL );
                                                // value = (double) get_volume_real_value( volume, voxx, voxy, voxz, 0, 0);
                                                fprintf(outputfile, "%lf\t", value );
                                        } else {
                                                printf("Voxel %d %d %d is outside %s\n", voxx, voxy, voxz, cur_minc);
                                                fprintf(outputfile, "\t");
                                        }
                        }
                        fprintf(outputfile, "\n");

                        delete_volume(volume);
                        
                } else {
                        keep_looping = 0;
                }
        }
        fclose(glimfile);
        fclose(outputfile);

        free( x );
        free( y );
        free( z );

        return(0);

}
