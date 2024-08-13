
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

import sys
import matplotlib         as mpl
import matplotlib.pyplot  as plt
import numpy              as np
import os
import argparse as ap
from readrawdata import readrawdata

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

ha2ev = 27.211396132 ## conversion hartree -> eV
bo2an = 0.52917721067



def option_parser():
    """
    Parse command line inputs for exciton weights

    Parse:
        directory
        files
        scale
        ymin
        ymax
        title
        grid
        nrows
        ncols

    :return input_options: Dictionary of parsed command line arguments
    """
    p = ap.ArgumentParser(description=\
                'Plot exciton weight data from single and multiple files as subplots.')

    help_directory = 'List of the directories in which the data to be plotted have to be found. If only one or no directory is specified, the data for the plots are taken from the same directory. Default value is the current directory.'

    help_files = 'List of file names containing data to plot. At least one file must be specified.'

    help_title = "Used as --title 'String as a title' assign a title to the plot."
    help_no_title = 'If present, it disables the writing of the title.'

    help_grid = 'If present, a grid is plotted in correspondence to the position of the major ticks.'

    help_scale = "One float corresponding to the scaling factor in the horizontal and vertical size of the plot appearence."

    help_ymin = "Minimum y-value to output from."
    help_ymax = "Maximum y-value to output from."

    help_nrows = "Number of rows in subplot."
    help_ncols = "Number of columns in subplot."
    help_sharex = "Whether the x-axes in subplot should share the same domain."
    help_sharey = "Whether the y-axes in subplot should share the same domain."
    #---------------------------------------------------------------------------

    p.add_argument('-d','--directory',
                   nargs = '*', default = ["."],
                   type = str, help = help_directory)

    p.add_argument('-f','--files',
                   nargs = '*', default = [],
                   type = str, help = help_files)

    p.add_argument('-y1','--ymin',
                   nargs = '*', default = [None, None],
                   type = float, help = help_ymin)

    p.add_argument('-y2','--ymax',
                   nargs = '*', default = [None, None],
                   type = float, help = help_ymax)

    p.add_argument('-s','--scale', default = 1.0,
                   type = float, help = help_scale)

    p.add_argument('-t','--title',
                   type = str, default = None, help = help_title)

    p.add_argument('-nt','--no_title', action='store_true', help = help_no_title)

    p.add_argument('-nr','--nrows', type=int, help = help_nrows)
    p.add_argument('-nc','--ncols', type=int,  help = help_ncols)


    # p.add_argument('-lp','--legend_position',
    #                type = str, help = help_legend_position,
    #                choices = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
    #                           'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center'],
    #                default = 'best')

    # p.add_argument('-nl','--no_legend', action='store_true', help = help_no_legend)

    p.add_argument('-g','--grid', action='store_true', help = help_grid)
    p.add_argument('-shx','--share_x', default=True, help = help_sharex)
    p.add_argument('-shy','--share_y', default=True, help = help_sharey)


    #---------------------------------------------------------------------------

    args = p.parse_args()
    input_options = {}

    input_options['directory'] = args.directory

    if ( len(args.files)==0 ):
        sys.exit("\n ERROR: At least a filename must be specified:\n\n"
                 +"        PLOT-exciton-weights.py -f [FILES [FILES ...]]\n")
    input_options['files'] = args.files

    input_options['ymin'] = None
    if ( len(args.yboundary) >= 1 ): input_options['ymin'] = args.yboundary[0]
    input_options['ymax'] = None
    if ( len(args.yboundary) >= 2 ): input_options['ymax'] = args.yboundary[1]

    input_options['scale'] = args.scale
    input_options['title'] = args.title
    input_options['no_title'] = args.no_title
    input_options['grid'] = args.grid
    input_options['nrows'] = args.nrows
    input_options['ncols'] = args.ncols
    input_options['share_x'] = args.share_x
    input_options['share_y'] = args.share_y

    return input_options




# narg = len(sys.argv)-1
# if( narg < 1):
#     print( "\n** ERROR: Must specify name of file and direction on command line.\n")
#     print( "** Usage: ", sys.argv[0], " <KS band structure width weights>")
#     sys.exit(0)
# if( not os.path.isfile( sys.argv[1])):
#     print( "\n** ERROR: Input file %s was not found." % sys.argv[1])
#     sys.exit(0)
# inname = sys.argv[1].strip()
# if( narg == 1):
#     scale = 1.0
# elif( narg == 2):
#     scale = float( sys.argv[2])
# elif( narg == 3):
#     ymin = float( sys.argv[2])
#     ymax = float( sys.argv[3])
#     scale = 1.0
# else:
#     ymin = float( sys.argv[2])
#     ymax = float( sys.argv[3])
#     scale = float( sys.argv[4])
# scale *= 2e2
# tmp = inname.split( '.')
# tmp[-1] = 'png'
# outname1 = inname+'.png'
# tmp[-1] = 'pdf'
# outname2 = inname+'.pdf'
# tmp[-1] = 'eps'
# outname3 = inname+'.eps'
# tmp[-1] = 'svg'
# outname4 = inname+'.svg'

# Remove or add image formats
output_formats: list = ["png", "eps", "pdf", "svg"]
def read_input(filename: str) -> tuple:
        """
        Reads KPATH and BANDLINE data.
        """
        bandwgt, bandwgtdim = readrawdata(filename)
        bandlin, bandlindim = readrawdata( '../BANDLINES.OUT')
        bandlines = bandlin[:,0,0]
        return bandwgt, bandwgtdim, bandlin, bandlindim, bandlines

def main(input_options):
    '''
    input:
    :input_options: dictionary that holds the input options parsed 
                    from the command line arguments
    '''
    directory = input_options['directory']
    files: list[str] = input_options['files']
    ymin: float = input_options['ymin'] 
    ymax: float = input_options['ymax'] 
    scale: float = input_options['scale'] 
    title: str = input_options['title']
    no_title: bool = input_options['no_title']
    grid: bool = input_options['grid']
    ncols: int = input_options['ncols']
    nrows: int = input_options['nrows']
    share_x: bool = input_options['share_x']
    share_y: bool = input_options['share_y']
    

    if (len(files) > nrows * ncols):
        print("WARNING: Number of plots greater than available subplot space. Please increase either nrows or ncols.")
        sys.exit(1)

    #########################
    # Settings for the plot #
    #########################

    figcolor = 'white'
    fig, axs = plt.subplots(nrows, ncols, sharex=share_x, sharey=share_y, figsize=(10,10))
    fig.patch.set_edgecolor(figcolor)
    fig.patch.set_facecolor(figcolor)

    mpl.rcParams['axes.linewidth'] = 3.0 # set the value globally
    mpl.rcParams['grid.linewidth'] = 1.5
    mpl.rcParams['xtick.labelsize'] = 30
    mpl.rcParams['ytick.labelsize'] = 30
    mpl.rcParams['axes.edgecolor'] = 'black'
    mpl.rcParams['axes.labelsize'] = '30'     # fontsize of the x any y labels
    mpl.rcParams['axes.labelcolor'] = 'black'
    mpl.rcParams['axes.axisbelow'] = 'True'   # whether axis gridlines and ticks are below
                                            # the axes elements (lines, text, etc)
    mpl.rcParams['legend.fontsize'] = 25
    plt.rcParams['xtick.major.pad'] = 10
    plt.rcParams['ytick.major.pad'] = 10

    plt.rcParams.update({'mathtext.default':'regular'})

    # Extract Data and Plot
    data = []
    num_plots = axs.flatten()
    counter = 0
    for row in axs:
        for col in row:
            col2 = col.twinx()
            output = read_input(files[counter])
            col.xaxis.grid( True, which='major', color='k', linestyle='-', linewidth=2)
            col.xaxis.set_label_position( 'bottom')
            col.set_xticks(output[4])
            col.set_xticklabels( ('W','L','$\Gamma$','X','W','K'))
            col.set_ylabel( 'Energy [eV]')
            for line in col.get_xticklines() + col.get_yticklines():
                line.set_markersize(10)
                line.set_markeredgewidth(2)
            xmin = np.amin( output[0][:,0,:])
            xmax = np.amax( output[0][:,0,:])
            col.set_xlim( [xmin, xmax])
            col.set_ylim( [ymin, ymax])

            # band-structure and weights
            for i in range( output[1][2]):
                col.plot( output[0][i,0,:], output[0][i,1,:], 'b', lw=3.0, zorder=10)
            for i in range( output[1][2]):
                col.scatter( output[0][i,0,:], output[0][i,1,:], s=(scale*output[0][i,2,:])**2, lw=3.0, edgecolor='r', facecolor='none', zorder=11)

            # Fermi level
            col.plot( [xmin, xmax], [0, 0], 'k', lw=3.0, ls='-')

            col2.set_ylim( col.get_ylim())
            col2.set_yticks( [0, 0])
            col2.set_yticklabels( ('$\\mathregular{E_{F}}$', ''))

            col.grid( True)
            plt.title( "", fontsize=mpl.rcParams['ytick.labelsize'], y=1.03)
            counter += 1
    #save file format
    for format in output_formats:
        fig.savefig(f"PLOT.{format}", format=format, dpi=300, bbox_inches='tight')
    # fig.savefig(outname3, format="eps",bbox_inches='tight')
    # fig.savefig( outname2, format='pdf', bbox_inches='tight')
    # fig.savefig(outname4, format="svg", bbox_inches='tight')
        # ax2 = ax1.twinx()

        #  axs[] = fig.add_axes( [0.17,0.1,0.75,0.8])
    ######################
    # Bandstructure plot #
    ######################

    
    # ax1.xaxis.grid( True, which='major', color='k', linestyle='-', linewidth=2)
    # ax1.xaxis.set_label_position( 'bottom')
   
    # ax1.set_xticklabels( ('W','L','$\Gamma$','X','W','K'))
    # ax1.set_ylabel( 'Energy [eV]')

    # Tick size
    
    # xmin = np.amin( bandwgt[:,0,:])
    # xmax = np.amax( bandwgt[:,0,:])
    # ax1.set_xlim( [xmin, xmax])
    # if( narg > 2): ax1.set_ylim( [ymin, ymax])

    # # band-structure and weights
    # for i in range( bandwgtdim[2]):
    #     ax1.plot( bandwgt[i,0,:], bandwgt[i,1,:], 'b', lw=3.0, zorder=10)
    # for i in range( bandwgtdim[2]):
    #     ax1.scatter( bandwgt[i,0,:], bandwgt[i,1,:], s=(scale*bandwgt[i,2,:])**2, lw=3.0, edgecolor='r', facecolor='none', zorder=11)

    # # Fermi level
    # ax1.plot( [xmin, xmax], [0, 0], 'k', lw=3.0, ls='-')
    # ax2.set_ylim( ax1.get_ylim())
    # ax2.set_yticks( [0, 0])
    # ax2.set_yticklabels( ('$\\mathregular{E_{F}}$', ''))

    # ax1.grid( True)
    # plt.title( "", fontsize=mpl.rcParams['ytick.labelsize'], y=1.03)

    # fig.savefig( outname1, format='png', dpi=300, bbox_inches='tight')
    # fig.savefig(outname3, format="eps",bbox_inches='tight')
    # fig.savefig( outname2, format='pdf', bbox_inches='tight')
    # fig.savefig(outname4, format="svg", bbox_inches='tight')
    #plt.show()
    sys.exit()
