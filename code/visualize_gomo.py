import pandas as pd
import argparse

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='Visualize GOMO motif enrichment file (XML)')   
    argparser.add_argument('--xml_from_gomo', type=str)
    args = argparser.parse_args()
    xml_from_gomo = args.xml_from_gomo
    gomo_plot_path = '/'.join(['/'.join(xml_from_gomo.split('/')[:-1]), 'gomo.png'])
    print(gomo_plot_path)
    
    df = pd.read_xml(xml_from_gomo)
    print(df.head())