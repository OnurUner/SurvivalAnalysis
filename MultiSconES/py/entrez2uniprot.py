# -*- coding: utf-8 -*-
import threading
from PyEntrezId import Conversion

class convert_thread(threading.Thread):
    def __init__(self, threadID, entrez_ids):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.entrez_ids = entrez_ids
        self.uniprot_ids = dict()       
        self.Id = Conversion('onurcan.uner@gmail.com')
        self.missing_ids = list()

    def run(self):
        print "Starting thread", self.threadID
        self.entrez2uniprot()
        print "Finished thread", self.threadID 
        
    def entrez2uniprot(self):
        for entrez_id in self.entrez_ids:        
            try:
                uniprot_id = self.Id.convert_entrez_to_uniprot(entrez_id)
                self.uniprot_ids[entrez_id] = uniprot_id
            except:
                self.missing_ids.append(entrez_id)
            
        