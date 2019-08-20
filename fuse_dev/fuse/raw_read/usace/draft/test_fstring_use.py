# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 14:06:32 2019

@author: jkinney
"""
my_etree_dict1 = {}
len_root_name_to_remove = len(xml_data.xml_tree.tag)
vertdatum = {'metadata/spref/vertdef/altsys/altdatum': 'altdatum'}
for key in iso_xml_path_to_baseattribute:
    if xml_data.xml_tree.findall(f'./{key[len_root_name_to_remove:]}'):
        if xml_data.xml_tree.findall(f'./{key[len_root_name_to_remove:]}') is None:
            my_etree_dict1[iso_xml_path_to_baseattribute[key]] = ''
        elif xml_data.xml_tree.findall(f'./{key[len_root_name_to_remove:]}') is list:  # check if list
            if len(xml_data.xml_tree.find(f'./{key[len_root_name_to_remove:]}')) > 0:
                my_etree_dict1[iso_xml_path_to_baseattribute[key]] = \
                    xml_data.xml_tree.find(f'./{key[len_root_name_to_remove:]}')[0].text
            else:
                my_etree_dict1[iso_xml_path_to_baseattribute[key]] = xml_data.xml_tree.find(
                    f'./{key[len_root_name_to_remove:]}').text
    else:
        my_etree_dict1[iso_xml_path_to_baseattribute[key]] = ''
for key in vertdatum:  #
    if xml_data.xml_tree.findall(f'./{key[len_root_name_to_remove:]}'):
        if isinstance(xml_data.xml_tree.findall(f'./{key[len_root_name_to_remove:]}'), list) == True:
        #if xml_data.xml_tree.findall(f'./{key[len_root_name_to_remove:]}') is list:  # check if list
            if len(xml_data.xml_tree.find(f'./{key[len_root_name_to_remove:]}')) > 0:
                if xml_data.xml_tree.find(
                        f'./{key[len_root_name_to_remove:]}') is None:  # Checks for NoneType object ('None')
                    my_etree_dict1['script: from_vert_key'] = ''
                    my_etree_dict1['from_vert_key'] = ''
                else:
                    my_etree_dict1[vertdatum[key]] = xml_data.xml_tree.find(
                        f'./{key[len_root_name_to_remove:]}').text
                    my_etree_dict1['from_vert_key'] = my_etree_dict1[vertdatum[key]]
            else:
                my_etree_dict1['from_vert_key'] = xml_data.xml_tree.find(f'./{key[len_root_name_to_remove:]}').text
        else:
            my_etree_dict1['from_vert_key'] = xml_data.xml_tree.find(f'./{key[len_root_name_to_remove:]}')
    else:
        my_etree_dict1['from_vert_key'] = ''
    my_etree_dict1['script: from_vert_key'] = my_etree_dict1['from_vert_key']
for x in xml_data.xml_tree.findall('.//ellips'):
    if xml_data.xml_tree.findall('.//ellips') is None:
        my_etree_dict1['ISO_ellips'] = ''
    else:
        my_etree_dict1['ISO_ellips'] = 'Exists'
        my_etree_dict1['ISO_xml'] = 'True'