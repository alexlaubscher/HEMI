# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 18:37:45 2018

@author: ePogue
"""

#!/usr/bin/env python
import threading, logging, time
import multiprocessing
import os
import sys
import time
import glob
from kafka import KafkaConsumer, KafkaProducer



#Create producer object

class Producer(threading.Thread):
    def __init__(self, msg):
        threading.Thread.__init__(self)
        self.msg=msg
        self.stop_event = threading.Event()
        
    def stop(self):
        self.stop_event.set()

    def run(self):
        #print(self.msg)
        producer = KafkaProducer(bootstrap_servers='localhost:9092')

        while not self.stop_event.is_set():
            producer.send('ppms', self.msg)

            time.sleep(1)

        producer.close()

class Consumer(multiprocessing.Process):
    def __init__(self):
        multiprocessing.Process.__init__(self)
        self.stop_event = multiprocessing.Event()
        
    def stop(self):
        self.stop_event.set()
       
    def run(self):
        consumer = KafkaConsumer(bootstrap_servers='localhost:9092',
                                 auto_offset_reset='earliest',
                                 consumer_timeout_ms=1000)
        consumer.subscribe(['my-topic'])

        while not self.stop_event.is_set():
            for message in consumer:
                print(message)
                if self.stop_event.is_set():
                    break

        consumer.close()
def findMostRecentFile():
    #find most recent file in log folder
    fName=max(glob.iglob('*.dat'),key=os.path.getmtime)
    return fName
        
def main():
    #find most recent file in log folder
    fname=findMostRecentFile()
    timeLastSaved=os.stat(fname).st_mtime
    
    #def kafkaExport()
        
    #Continuously check if most recent file and current file (fname) have changed.
    #If current file has changed, send to Kafka
    
    #Initially post the last line
    running=True
    timeLastSavedCurrent=os.stat(fname).st_mtime
    time.sleep(1)  
    try:
        f=open(fname)
        lineList=f.readlines()
        f.close()
        lastLine=lineList[len(lineList)-1]
        #Parse the last line for kafka
        p=Producer(lastLine)
        p.start()
        time.sleep(10)
        p.stop()
        p.join()
    finally:
        f.close()
        
        
    while running==True:
        #check if file has changed
        timeLastSavedCurrent=os.stat(fname).st_mtime
        time.sleep(1)
        if timeLastSaved != timeLastSavedCurrent:
            timeLastSaved=timeLastSavedCurrent
            try:
                f=open(fname)
                lineList=f.readlines()
                f.close()
                lastLine=lineList[len(lineList)-1]
                #Parse the last line for kafka
                p.msg=lastLine
                p.start()
                time.sleep(10)
                p.stop()
                p.join()
            finally:
                f.close()
    
                
        #check if more recent file created, updating accordingly
        latestFile=findMostRecentFile()
        if latestFile!=fname:
            fname=latestFile
            timeLastSaved=os.stat(fname).st_mtime
        #Allow it to be read
        #c=Consumer()
        #c.start()
        #time.sleep(10)
        #c.stop()
        #c.join()
        
if __name__ == "__main__":
    logging.basicConfig(
        format='%(asctime)s.%(msecs)s:%(name)s:%(thread)d:%(levelname)s:%(process)d:%(message)s',
        level=logging.INFO
        )
    main()