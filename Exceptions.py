class HttpForbidden(RuntimeError):
  pass



class HttpGatewayTimeout(RuntimeError):
  pass



class HttpInternalServerError(RuntimeError):
  pass



class RemoteApiRequestsExceeded(RuntimeError):
  pass
