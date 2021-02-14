from pymongo import MongoClient
def mongoConnect():
    """
    Connect to eggNOG database using pymongo
    """
    client = MongoClient('10.0.3.1', 27017, maxPoolSize=10)
    db = client.freeze11
    return client, db
