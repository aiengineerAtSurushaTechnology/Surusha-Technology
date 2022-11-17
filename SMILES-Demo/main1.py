# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 09:35:51 2022

@author: Acer
"""

#%%
import streamlit as st
import requests

from last import run_search

def main():
    menu = ['Home','Search','EDA','ML/DL','About']
    choice = st.sidebar.selectbox("Menu",menu)
    
    st.title("Chemical Reaction Prediction")
    
    if choice == 'Home':
        st.subheader('Home')
        st.success(""":pray: This is the Demo/prototype ML web app which  is used to present a baseline layout of our final 
                   app that is how minimum features our web app will have and what it can do. In this demo web app we 
                   separate the web app into 5 parts they can be see by sidebar, every option on the sidebar have its specific 
                   name which represnt what work that secssion will do.:clap:""")
        
    elif choice == 'Search':
        st.subheader("Sreach for SMILES 3D-Graph")
        st.subheader("To print the 3D-molecular represntation we use RDkit and Py3Mol libraries.")
        run_search()
    elif choice == "EDA":
        st.subheader("The Exploratory Data Analysis will representated in this section of the web app")
        pass
    elif choice == "ML/DL":
        st.subheader("The ML or DL model will deploy in this area of the web app to predic the chemical reactions.")
        pass
    else:
        st.subheader("About")
        st.success("""This Silnelia web app is created by Rohit Kumar :snowman:. Who is presently a Junior Data Scientist at 
                 Surusha Technology Ltd.:sunny:. This is the demo work what we will try to do in future this demo have limited functionality 
                 and work exposure. To creat this demo I used python's streamlit library which is very good for 
                 this kind of demo work because we can create this type of web app very easily with the help of this library.
                 In this demo work we first creat a home page after that we create a EDA page which may be in future change as 
                 to give the  SMILES structures of the molecules for the chemical formula or just vice-versa.
                 In the next page we do our ML work that is our ML algprithms are deploy in this page and the user may 
                 use them for there work. At last it is the about page.:snowflake:.""")
        
if __name__ == '__main__':
    main()