
#include <mrc_json.h>
#include <json-builder.h>

#include <string.h>
#include <stdio.h>
#include <assert.h>

int
main(int argc, char **argv)
{
  json_settings settings = {};
  settings.value_extra = json_builder_extra;  /* space for json-builder state */

  char json_str[] = "{ \"i\": 2.0, \"i3\": [ 2, 3, 4] }";

  char error[json_error_max];
  mrc_json_t json = {
    .value = json_parse_ex(&settings, json_str, strlen(json_str), error),
  };
  assert(json.value);

  char *buf = mrc_json_to_string(json);
  printf("%s\n\n", buf);
  free(buf);
  
  mrc_json_print(json, 0);

  return 0;
}
