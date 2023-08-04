class HttpForbidden(RuntimeError):
  pass



class HttpRequestTimeout(RuntimeError):
  pass



class HttpInternalServerError(RuntimeError):
  pass



class RemoteApiRequestsExceeded(RuntimeError):
  pass
