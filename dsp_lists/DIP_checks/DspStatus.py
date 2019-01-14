from requests import Session
from requests.auth import HTTPBasicAuth  # or HTTPDigestAuth, or OAuth1, etc.
from zeep import Client
from zeep.transports import Transport

session = Session()
session.auth = HTTPBasicAuth('pporras', 'yourpassword')

client = Client('https://imexcentral.org/icentral/ws-v20?wsdl',
                transport=Transport(session=session))

identifier_type = client.get_type('ns0:identifier')
publication_type = client.get_type('ns0:publication')

with client.settings():
    with open('dsp.txt') as pmids:
        with open("DspStatus" + ".txt", "w") as f:
            f.write("Pubmed\tOwner\tStatus\tAdmin Users\tAdmin Groups\t\n")
            line = pmids.readline()
            cnt = 1
            while line:
                ac = line.strip()

                identifier = identifier_type(ns='pmid', ac=ac)
                response = client.service.getPublicationById(identifier=identifier)

                users = ''
                groups = ''

                if response.adminUserList:
                    users = ','.join(response.adminUserList.user)

                if response.adminGroupList:
                    groups = ','.join(response.adminGroupList.group)

                f.write("{}\t{}\t{}\t{}\t{}\n".format(ac, response.owner, response.status, users, groups))
                line = pmids.readline()
                cnt += 1
