import requests
import sys

headers = {"api-key":"ff352a2affacf64753689dd603b5b44e"}

def get(path, params=None):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically
    return r


snapshot = sys.argv[1]
subfind = sys.argv[2]

url = f'https://www.tng-project.org/api/TNG300-1/snapshots/{snapshot}/subhalos/{subfind}/'

print(url)


